from assembler.gaf_reader import GafReader
#from assembler.builder import AnchorDictionary
from assembler.aligner import AlignAnchor
import assembler.aligner as aligner_module
import assembler.parser as parser
from assembler.helpers import traced_functions
import time
import multiprocessing
from sys import stderr
import os
from assembler.config import settings
from collections import defaultdict
from assembler.read import Read
if not os.environ.get("MEMORY_PROFILE"):
    try:
        from line_profiler import profile
    except ImportError:
        def profile(func):
            return func

# Ensure we use fork mode for true copy-on-write behavior
# (on Linux, this is the default, but we make it explicit for clarity)
try:
    multiprocessing.set_start_method('fork', force=False)
except RuntimeError:
    # Already set, which is fine
    pass

# Global objects to hold shared data for worker processes
# Note: Set these before creating the multiprocessing pool to leverage fork()'s copy-on-write
shared_align_anchor = None
shared_gaf_chunks = None

def nested_dd_factory():
    # Deliberately undecorated: this is used as a defaultdict factory, and a decorated
    # (closure-wrapped) factory gets embedded in the dict that process_gaf_chunk returns —
    # multiprocessing.Pool then fails to pickle that dict on the way back to the parent.
    return defaultdict(list)

@profile
def init_worker():
    """
    Initializer for each worker process in the pool.
    On Linux with fork, the global shared_align_anchor is already available via copy-on-write.
    """
    # No need to do anything - fork() gives us access to the parent's memory
    pass


def process_gaf_chunk(chunk_idx: int) -> dict:
    """
    Picklable top-level dispatcher for a chunk of GAF lines.
    Kept undecorated so multiprocessing.Pool can pickle it under every profiling mode —
    memory_profiler wraps @profile-decorated functions in a closure that pickle can't
    serialize by reference, which used to force MEMORY_PROFILE runs to skip Pool entirely.

    Under MEMORY_PROFILE, `python -m memory_profiler` reports results exactly once, in a
    finally block wrapping the whole script, in whichever process reaches it — that's the
    parent, never a forked worker, since workers just get reaped when the Pool closes. So
    under MEMORY_PROFILE this builds a fresh, worker-local LineProfiler and writes its own
    report per worker instead of relying on the module-level @profile / global profiler.
    Under every other mode, _process_gaf_chunk_impl already carries the real @profile
    decorator (bound at import time) and this just calls it directly.

    Takes an integer index into shared_gaf_chunks (set as a module global before the Pool
    is created, so workers inherit the actual chunk data via fork copy-on-write) rather than
    the chunk's lines directly — pool.map otherwise pickles whatever it's given through the
    IPC pipe to send to each worker, and the chunk lines are the largest single payload in
    this pipeline (previously ~575MB of pickling across all workers combined).
    """
    gaf_chunk_lines = shared_gaf_chunks[chunk_idx]
    if os.environ.get("MEMORY_PROFILE"):
        import memory_profiler
        prof = memory_profiler.LineProfiler(backend=memory_profiler.choose_backend('psutil'))
        # Nested calls a fresh LineProfiler wouldn't otherwise see: processGafLine walks
        # each node, calling verify_path_concordance then verify_sequence_agreement per
        # anchor hit, which itself calls parse_cs_tag's caller upstream in parser.py.
        targets = [
            (parser, "processGafLine"),
            (parser, "parse_cs_tag"),
            (aligner_module, "verify_sequence_agreement"),
            (aligner_module, "verify_path_concordance"),
        ]
        with traced_functions(prof, targets):
            result = prof(_process_gaf_chunk_impl)(gaf_chunk_lines)
        with open(f"worker_gaf_mem_{os.getpid()}.txt", "w") as report_f:
            memory_profiler.show_results(prof, stream=report_f)
        return result
    return _process_gaf_chunk_impl(gaf_chunk_lines)


def _process_gaf_chunk_impl(gaf_chunk_lines: list[str]) -> dict:
    """
    Worker function to process a chunk of GAF lines.
    This function is executed in a separate process.
    """
    if hasattr(profile, '_profile') and profile._profile is not None:
        profile.enable()
    global shared_align_anchor

    # Initialize local dictionaries to store results for this chunk.
    local_anchor_reads_dict = defaultdict(nested_dd_factory)
    local_bp_matched_reads = defaultdict(list)
    local_path_matched_reads = defaultdict(list)
    local_reads_processed_dict = {} # {read_name: processed_line_data}

    # Worker-local cache of (node_handle, length) per node ID. Reads within a chunk share
    # most of their path, so the vast majority of node visits are repeats — this trades
    # 3 C-extension calls (has_node/get_handle/get_length) for 1 dict lookup on a hit.
    node_cache = {}

    t0 = time.time()
    # Process each line in the assigned GAF chunk.
    for line in gaf_chunk_lines:
        processed_line_data = parser.processGafLine(line)

        if shared_align_anchor.read_id_map:
            if not processed_line_data:
                continue
            read_name = processed_line_data[0]
            read_id = shared_align_anchor.read_id_map.get(read_name)
            if read_id is None:
                continue # Skip GAF entries for reads not in the FASTA file
            processed_line_data[0] = read_id

        if processed_line_data:
            if settings.OUTPUT_LOGGING_FILES:
                local_reads_processed_dict[processed_line_data[0]] = processed_line_data
            # remove mapq (index MAP_Q_ID) and div (index DIV_ID) from processed_line_data
            processed_line_data = processed_line_data[:4] + processed_line_data[6:]
        
            # Call the refactored processGafLine on the shared object
            # This is a read-only operation on shared_align_anchor
            result, current_read = shared_align_anchor.processGafLine(processed_line_data, node_cache=node_cache)
                        
            for (sentinel, i), reads in result["anchor_reads"].items():
                local_anchor_reads_dict[sentinel][i].extend(reads)
            
            for (sentinel, i), reads in result["bp_matched_reads"].items():
                local_bp_matched_reads[(sentinel, i)].extend(reads)

            if settings.OUTPUT_LOGGING_FILES:
                for (sentinel, i), reads in result["path_matched_reads"].items():
                    local_path_matched_reads[(sentinel, i)].extend(reads)

    if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
        print(f" ..Processed {len(gaf_chunk_lines)} lines in {time.time()-t0:.2f}s", file=stderr)

    if hasattr(profile, '_profile') and profile._profile is not None:
        profile._profile.dump_stats(f"worker_gaf_{os.getpid()}.lprof")
    # Return the collected results from this worker.
    return {
        "anchor_reads_dict": local_anchor_reads_dict,
        "bp_matched_reads": local_bp_matched_reads,
        "path_matched_reads": local_path_matched_reads,
        "reads_processed": local_reads_processed_dict
    }


if not os.environ.get("MEMORY_PROFILE"):
    # line_profiler mode: decorate statically, same as every other @profile use in this file.
    # Left undecorated under MEMORY_PROFILE — process_gaf_chunk wraps the raw function itself
    # with a fresh, worker-local LineProfiler; wrapping an already-decorated function here
    # would make that fresh profiler trace the wrapper's code instead of the real one.
    _process_gaf_chunk_impl = profile(_process_gaf_chunk_impl)


class Orchestrator:

    @profile
    def __init__(
        self, dictionary_path: str, graph_path: str, gaf_path: str, fasta_path: str, threads: int, read_id_map: dict = None
    ):
        """
        It initiailzes the AlignAnchor object with the packedgraph path and the dictionary generated by the assembler.builder.AnchorDictionrary object.
        It initializes the GafReader object that reads the gaf file.

        Parameters
        ----------
        sentinel_to_anchor_dictionary: dictionary
            the dctionary associating sentinels and anchors
        graph_path: string
            The filepath of the packedGraph object
        gaf_path:
            The filepath of the gaf alignment file
        fasta_path: string
            The filepath of the reads fasta file
        """
        self.align_anchor = AlignAnchor(threads=int(threads), read_id_map=read_id_map)
        t0 = time.time()
        self.align_anchor.build(dictionary_path, graph_path)    # graph is loaded here once!
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"AlignAnchor built in {time.time()-t0:.2f}s", file=stderr)
        self.align_anchor.readFasta(fasta_path)
        self.gaf_path = gaf_path
        self.threads = int(threads)
        self.total_reads_in_gaf = 0

    @profile
    def _chunk_gaf_file(self, gaf_path: str, num_chunks: int) -> list:
        """
        Reads a GAF file and splits its lines into a specified number of chunks for parallel processing.
        """
        with open(gaf_path, "r") as f:
            lines = f.readlines()
        
        self.total_reads_in_gaf = len(lines)
        if not lines:
            return []

        chunk_size = (len(lines) + num_chunks - 1) // num_chunks
        return [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)] # output: [[line1, line2, ...], [line6, line7, ...], ...]

    @profile
    def process(self, out_prefix: str, debug_file=None):
        """
        Orchestrates the processing of the GAF file, either in a single thread or in parallel.
        """
        t0 = time.time()
        
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(f"Processing GAF file in parallel with {self.threads} threads...", file=stderr)
        
        # Set globals before forking to leverage copy-on-write (avoids pickling)
        global shared_align_anchor, shared_gaf_chunks
        shared_align_anchor = self.align_anchor

        os.environ["LINE_PROFILE"] = "1"
        # Divide the GAF file into chunks
        shared_gaf_chunks = self._chunk_gaf_file(self.gaf_path, self.threads)

        # process_gaf_chunk is an undecorated dispatcher (see its docstring), so it pickles
        # fine under Pool.map regardless of profiling mode, including MEMORY_PROFILE. It only
        # takes an integer chunk index now — the actual chunk lines are inherited by each
        # worker via fork COW from shared_gaf_chunks, not pickled through the IPC pipe.
        with multiprocessing.Pool(processes=self.threads, initializer=init_worker) as pool:
            results = pool.map(process_gaf_chunk, range(len(shared_gaf_chunks)))
        
        if settings.DEBUG:
            print("Merging results from worker processes...", file=stderr)
        # Prepare reads_processed TSV: remove old file once before appending
        reads_processed_path = f"{out_prefix}.reads_processed.tsv" if settings.OUTPUT_LOGGING_FILES else None
        if reads_processed_path and os.path.exists(reads_processed_path):
            os.remove(reads_processed_path)
        for result_dict in results:
            self.align_anchor.merge_results(result_dict, reads_processed_path)
        
        total_time_for_gaf_processing = time.time() - t0
        
        if settings.DEBUG or settings.PRINT_RUNTIME_LOGS:
            print(
                f"GAF processing finished in {total_time_for_gaf_processing:.2f}s with {self.threads} threads",
                file=stderr,
            )

        # Run the dump_valid_anchors method which runs the unreliable snarl filtering and the anchor extensions
        
        kwargs = {
            "extended_out_file_path": f"{out_prefix}.extended.jsonl",
            "reliable_snarls_out_file_path": f"{out_prefix}.reliable_snarls.tsv",
            "pre_reliable_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.pre_reliable.tsv",
        }

        if settings.OUTPUT_LOGGING_FILES:
            kwargs.update({
                "path_matched_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.path_matched.tsv",
                "seq_matched_sizes_out_file_path": f"{out_prefix}.subgraph.sizes.seq_matched.tsv",
                "anchor_read_tracking_file_path": f"{out_prefix}.read_drop_tracking.jsonl",
                "independent_anchor_read_tracking_file_path": f"{out_prefix}.independent_ext_tracking.jsonl",
                "snarl_variant_type_out_file_path": f"{out_prefix}.snarl_variant_type.jsonl",
                "snarl_compatibility_out_file_path": f"{out_prefix}.snarl_compatibility.jsonl",
                "snarl_common_reads_out_file_path": f"{out_prefix}.snarl_2_snarl_common_reads.jsonl",
                "snarl_read_partitions_out_file_path": f"{out_prefix}.snarl_2_snarl_read_partitions.jsonl",
                "snarl_coverage_out_file_path": f"{out_prefix}.snarl_coverage.jsonl",
                "snarl_allelic_coverage_out_file_path": f"{out_prefix}.snarl_allelic_coverage.jsonl",
                "snarl_coverage_extended_out_file_path": f"{out_prefix}.snarl_coverage_extended.jsonl",
                "snarl_allelic_coverage_extended_out_file_path": f"{out_prefix}.snarl_allelic_coverage_extended.jsonl",
                "binomial_pairs_out_file_path": f"{out_prefix}.binomial_pairs.tsv"
            })
            self.align_anchor.dump_valid_anchors(**kwargs)
            self.align_anchor.dump_snarls_and_anchors_in_reads_dict(f"{out_prefix}.snarls_and_anchors_in_reads.jsonl")
        
        else:
            # Just dump the extended valid anchors JSON
            self.align_anchor.dump_valid_anchors(**kwargs)

        # Always emit the extended subgraph size TSV, independent of OUTPUT_LOGGING_FILES.
        # This is a lightweight summary artifact that downstream steps may rely on.
        out_file = f"{out_prefix}.subgraph.sizes.extended.tsv"
        out_dir = os.path.dirname(out_file)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        self.dump_dict_size_extended(out_file)


    @profile
    def dump_dictionary_with_counts(self, out_file: str):
        """
        It dumps the positioned anchor dictionary by json
        """
        self.align_anchor.dump_dictionary_with_reads_counts(out_file)

    @profile
    def dump_dict_size_extended(self, out_file: str):
        """
        It dumps the anchors by json
        """
        self.align_anchor.print_extended_anchor_info(out_file) 

    @profile
    def dump_bandage_csv_extended(self, out_file: str):
        """
        It dumps CSV with node and colour of all anchor nodes
        """
        self.align_anchor.print_sentinels_for_bandage(out_file) 
