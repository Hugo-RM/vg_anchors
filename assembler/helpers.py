from sys import argv, stderr, exit

import os
if not os.environ.get("MEMORY_PROFILE"):
    try:
        from line_profiler import profile
    except ImportError:
        def profile(func):
            return func
import json
from collections import defaultdict
from assembler.anchor import Anchor
from assembler.config import settings
import gzip
from contextlib import contextmanager
from Bio import SeqIO

@profile
def reverse_complement(string) -> str:
    rev_str = string[::-1]
    r_c = ""
    for el in rev_str:
        if el == "A":
            r_c += "T"
        if el == "C":
            r_c += "G"
        if el == "G":
            r_c += "C"
        if el == "T":
            r_c += "A"
    return r_c

@contextmanager
def open_fastq(filename):
    try:
        if filename.endswith(".gz"):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'r')
        try:
            yield f
        finally:
            f.close()
    except IOError as e:
        if settings.DEBUG:
            print(f"Error opening file {filename}: {e}")
        raise

@profile
def fastq_lines(in_fastqs):
    for fname in in_fastqs:
        if settings.DEBUG:
            print(fname,flush=True)
        with open_fastq(fname) as f:
            yield from f

@profile
def fastq_entries(fastq_lines_iter):
    """Generator that yields complete FASTQ entries"""
    while True:
        try:
            header = next(fastq_lines_iter)
            sequence = next(fastq_lines_iter)
            plus_line = next(fastq_lines_iter)
            quality = next(fastq_lines_iter)
            
            yield {
                'header': header.strip().split('\t')[0][1:],
                'sequence': sequence.strip(),
                'plus_line': plus_line.strip(),
                'quality': quality.strip()
            }
        
        except StopIteration:
            break


# Function to get complement
@profile
def complement(seq):
    # Define complement dictionary
    complement_map = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement_map)

# Function to get reverse complement
@profile
def rev_c(seq):
    return complement(seq)[::-1]

@profile
def extract_sequence(fasta_file, read_id):
    """
    get sequence for read from fasta
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == read_id:
            return str(record.seq)
    return None  # Return None if read_id is not found


@contextmanager
def traced_functions(prof, targets):
    """
    Temporarily wrap each (owner, attr_name) target with `prof` (a memory_profiler
    LineProfiler) so they all report into the same profile, then restore the originals.
    owner is a module or a class; attr_name is looked up on it with getattr/setattr, which
    works for plain module-level functions and for class methods (accessed unbound via the
    class, Python's descriptor protocol re-binds them to `self` normally when called).

    Used by the GAF/snarl worker dispatchers under MEMORY_PROFILE: a fresh LineProfiler
    only tracks the exact function object it's given, not anything that function calls, so
    without this, only the dispatcher's own top-level function shows up per worker.
    """
    originals = [(owner, name, getattr(owner, name)) for owner, name in targets]
    try:
        for owner, name, original in originals:
            setattr(owner, name, prof(original))
        yield
    finally:
        for owner, name, original in originals:
            setattr(owner, name, original)

if __name__ == "__main__":
    # verify_anchors_validity(argv[1], argv[2], argv[3])
    #anchors_shasta = argv[1]
    anchors_pos_dict = argv[1]
    out_png = argv[2]
    title = argv[3]

    plot_count_histogram(anchors_pos_dict, out_png + "count.png")

    plot_anchor_count_genome_distribution(anchors_pos_dict, out_png + "position_count.png", title
    )
