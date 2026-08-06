from sys import stderr
import re
from array import array
from itertools import accumulate
from assembler.config import settings

import os
if not os.environ.get("MEMORY_PROFILE"):
    try:
        from line_profiler import profile
    except ImportError:
        def profile(func):
            return func

# Precompiled at import time — reused for every call
_CS_OPS = re.compile(r'[+\-=][ACGTNacgtn]+|:[0-9]+|\*[ACGTNacgtn]{2}')

"""
This functions process the gaf alignment file from Giraffe HiFi and return the necessary data
structured for "fast" anchor queries.
- process_line takes start/end of the alignment in the path, the path and the cs:Z optional tag.
It parses the path field to record nodes and orientations, stored in two numpy arrays (int and bool)

- parse_cs_line parses the cs:Z field to structure the cigar information to verify the basepair alignment between the path sequence and the read 

For gaf tags see: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10
"""


def processGafLine(gaf_line: str):
    """
    It parses a GAF line to extract and structure useful tags and returns them in a list

    Parameters
    ----------
    gaf_line : string
        a gaf line stripped of whitespaces at the beginning and end

    Returns
    -------
    list
        list of processed tags:
        read_name : string
        read_len : int
        read_start : int - start of the alignment on the read
        relative_strand : bool (True if +, False else)
        path_start : int - start of the alignment in the path sequence
        path_end : int - end of the alignment in the path sequence
        nodes_list : list - of node_ids of the nodes walked by the path
        orientation_list : list - of node orientations of the nodes walked by the path
        cs_line : list - succession of tuples describing the cigar
        cum_path : array.array('q') - prefix sums of path deltas, length len(cs_line)+1, starting with 0
        cum_seq : array.array('q') - prefix sums of seq deltas, length len(cs_line)+1, starting with 0
    """

    line_elements = gaf_line.split()
    #print(line_elements[1:2])
    if not(line_elements[1].isnumeric()):
        #print('Not numeric')
        del line_elements[1:3]
    # First verify that the gaf line contains an usable alignment
    if (len(line_elements) in [settings.EXPECTED_GAF_TAGS, settings.EXPECTED_GAF_TAGS - 1]) and int(line_elements[settings.MAP_Q_ID]) >= settings.EXPECTED_MAP_Q:
        # extract needed tags
        read_name = line_elements[settings.READ_NAME_ID]
        read_len = int(line_elements[settings.READ_LEN])
        read_start = int(line_elements[settings.READ_START_ID])  # Will be used to initialize `walked_length` in the aligner.processGafLine function
        mapq = int(line_elements[settings.MAP_Q_ID])
        #relative_strand = True if line_elements[RELATIVE_STRAND_ID] == "+" else False
        path_start = int(line_elements[settings.PATH_START_ID])
        path_end = int(line_elements[settings.PATH_END_ID])
        # if the div tag is present, extract it, otherwise set it to 0
        div = 0
        if len(line_elements) == settings.EXPECTED_GAF_TAGS:
            if line_elements[settings.DIV_ID].startswith("dv:f:"):
                div = float(line_elements[settings.DIV_ID].split("dv:f:")[1])
        
        # decompose the path into nodes and orientations arrays
        nodes_list = []
        orientation_list = []

        curr_node_string = ""
        for char in line_elements[settings.PATH_ID]:
            if char in "><":
                orientation_list.append(True if char == ">" else False)
                if len(curr_node_string) != 0:
                    nodes_list.append(int(curr_node_string))
                    curr_node_string = ""
            else:
                curr_node_string += char
        if len(curr_node_string) != 0:
            nodes_list.append(int(curr_node_string))
        
        # FIXME: We don't need to output relative strand here. We are calculting that in the verify_path_concordance function.
        count_positive_orientation_nodes = orientation_list.count(True)
        relative_strand = True if count_positive_orientation_nodes > (len(orientation_list) / 2) else False

        # decompose the cs tag into alignment steps + precomputed cumulative offsets
        if len(line_elements[settings.CS_TAG_ID]) > settings.MIN_CS_LEN:
            cs_line, cum_path, cum_seq = parse_cs_tag(line_elements[settings.CS_TAG_ID])
        else:
            print("ERROR IN CS LINE.",flush=True, file=stderr)
            return None

        return [
            read_name,
            read_len,
            read_start,
            relative_strand,
            mapq,
            div,
            path_start,
            path_end,
            nodes_list,
            orientation_list,
            cs_line,
            cum_path,   # after handler slice: index 9 = CUM_PATH_POSITION
            cum_seq,    # after handler slice: index 10 = CUM_SEQ_POSITION
        ]

    if settings.DEBUG:
        print(f"ERROR: {len(line_elements)} =? {settings.EXPECTED_GAF_TAGS} _ {int(line_elements[settings.MAP_Q_ID])} =? {settings.EXPECTED_MAP_Q}",flush=True, file=stderr)
    return None


def parse_cs_tag(cs_string: str):
    """
    Parses the cs tag string into a list of alignment steps and cumulative position offsets.
    Each step is a tuple of (flag, val) where flag is the operation character and val is the length in bp.
    For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10

    Parameters
    ----------
    cs_string : string
        a string spelling the cs tag in the gaf

    Returns
    -------
    ops : list of (flag, val) tuples
        list of operations and bp movement
    cum_path : array.array('q')
        prefix sums of path deltas, length len(ops)+1, starting with 0.
        array.array instead of list: these are built once and only ever
        read via bisect/indexing afterward (never mutated), and every
        element is a plain int, so the per-element PyLong object overhead
        of a list buys nothing here — array.array stores them as packed
        8-byte C longs instead, cutting memory ~4x for these hot,
        per-read structures.
    cum_seq : array.array('q')
        prefix sums of seq deltas, length len(ops)+1, starting with 0. Same
        array.array rationale as cum_path.
    """
    # flag characters used to represent the basepair alignment
    # = : identical sequence, spelled [ACGTN]+
    # + : insertion to the reference, spelled [ACGTN]+
    # - : deletion to the reference, spelled [ACGTN]+
    # : : identical sequence, length [0-9]+
    # * : substitution (reference to query) [acgtn][acgtn]
    # example: cs:Z::6724+T:581+A:1027-G:2962
    ops = []
    path_deltas = []
    seq_deltas = []
    # Bound methods hoisted to locals: avoids repeated attribute lookup + method
    # dispatch on ops/path_deltas/seq_deltas for every cs operation in the string.
    append_op = ops.append
    append_path_delta = path_deltas.append
    append_seq_delta = seq_deltas.append
    for m in _CS_OPS.finditer(cs_string):
        op = m.group()
        flag = op[0]
        if flag == ':':
            val = int(op[1:])
            append_op((':', val)); append_path_delta(val); append_seq_delta(val)
        elif flag == '*':
            append_op(('*', 1)); append_path_delta(1); append_seq_delta(1)
        elif flag == '+':
            val = len(op) - 1
            append_op(('+', val)); append_path_delta(0); append_seq_delta(val)
        elif flag == '-':
            val = len(op) - 1
            append_op(('-', val)); append_path_delta(val); append_seq_delta(0)
        else:  # '='
            val = len(op) - 1
            append_op(('=', val)); append_path_delta(val); append_seq_delta(val)
    return ops, array('q', accumulate(path_deltas, initial=0)), array('q', accumulate(seq_deltas, initial=0))


if not os.environ.get("MEMORY_PROFILE"):
    # Left undecorated above and applied conditionally here: under MEMORY_PROFILE, worker
    # dispatchers (handler.py, aligner.py) wrap these with their own fresh, per-worker
    # LineProfiler via helpers.traced_functions(). A static @profile here would resolve to
    # the CLI's single global profiler instead, and the dispatcher would end up wrapping
    # that wrapper — tracking memory_profiler's own internals instead of this file's code.
    processGafLine = profile(processGafLine)
    parse_cs_tag = profile(parse_cs_tag)
