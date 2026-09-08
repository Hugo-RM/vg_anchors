from sys import stderr
import re
from itertools import accumulate
from assembler.config import settings

"""
This functions process the gaf alignment file from Giraffe HiFi and return the necessary data
structured for "fast" anchor queries.
- process_line takes start/end of the alignment in the path, the path and the cs:Z optional tag.
It parses the path field to record nodes and orientations, stored in two numpy arrays (int and bool)

- parse_cs_line parses the cs:Z field to structure the cigar information to verify the basepair alignment between the path sequence and the read 

For gaf tags see: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10
"""

# Precompiled at import time — reused for every call
_CS_OPS = re.compile(r'[+\-=][ACGTNacgtn]+|:[0-9]+|\*[ACGTNacgtn]{2}')


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
        cum_path : list - prefix sums of path deltas, length len(cs_line)+1, starting with 0
        cum_seq : list - prefix sums of seq deltas, length len(cs_line)+1, starting with 0
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
    Parses the cs tag string into a list of alignment steps, plus the cumulative path and
    sequence offsets after each step. Each step is a tuple of (flag, length) where flag is
    the operation character and length is the size of the operation in basepairs.
    For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10

    cum_path and cum_seq let verify_sequence_agreement binary search directly to the cs
    step an anchor starts in, instead of walking every step of the cs tag from the
    beginning for every anchor it checks against this read.

    Parameters
    ----------
    cs_string : string
        a string spelling the cs tag in the gaf

    Returns
    -------
    ops : list of (flag, length) tuples
        the cs tag operations, in order
    cum_path : list[int]
        running path position after each step (prefix sums of path deltas), length
        len(ops)+1, starting with 0. Only ever increases, so it can be binary searched.
    cum_seq : list[int]
        running read-sequence position after each step (prefix sums of seq deltas), same
        length and starting value as cum_path.
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
    for match in _CS_OPS.finditer(cs_string):
        op_str = match.group()
        flag = op_str[0]
        if flag == ':':
            length = int(op_str[1:])
            ops.append((':', length)); path_deltas.append(length); seq_deltas.append(length)
        elif flag == '*':
            ops.append(('*', 1)); path_deltas.append(1); seq_deltas.append(1)
        elif flag == '+':
            length = len(op_str) - 1
            ops.append(('+', length)); path_deltas.append(0); seq_deltas.append(length)
        elif flag == '-':
            length = len(op_str) - 1
            ops.append(('-', length)); path_deltas.append(length); seq_deltas.append(0)
        else:  # '='
            length = len(op_str) - 1
            ops.append(('=', length)); path_deltas.append(length); seq_deltas.append(length)
    return ops, list(accumulate(path_deltas, initial=0)), list(accumulate(seq_deltas, initial=0))
