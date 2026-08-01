from sys import stderr
from itertools import accumulate

import os
if not os.environ.get("MEMORY_PROFILE"):
    try:
        from line_profiler import profile
    except ImportError:
        def profile(func):
            return func

class Anchor:

    @profile
    def __init__(self) -> None:
        self._nodes: list = []
        self._node_ids = None  # lazy-build sentinel for build_lookup_cache(); see verify_path_concordance
        self.snarl_id: int = 0
        self.genomic_position: int = 0
        # self.baseparilength: int = 0
        self.basepairlength: int = 0
        self.sentinel_length: int = 0
        self.num_sequences: int = 0
        self.chromosome: str = ""
        self.reference_paths_covered: list = []
        # self.path_matched_reads: list = []
        self.bp_matched_reads: list = []        # [READ_ID, READ_STRAND, ANCHOR_START, ANCHOR_END, MATCH_LIMIT, CS_LEFT_AVAIL, CS_RIGHT_AVAIL]
        self._reads: list = []
        self.bp_occupied_start_node = 0       # basepairs occupied by the leftmost node in the anchor (this is independent of anchor orientation, which means that node with lowest node_id is considered leftmost)
        self.bp_occupied_end_node = 0
        self.path_orientation: bool = None

    @profile
    def copy_from_anchor(self, other_anchor):
        for attr_name, attr_value in other_anchor.__dict__.items():
            setattr(self, attr_name, attr_value)

    @profile
    def add_reference_path(self, path):
        self.reference_paths_covered.append(path)

    @profile
    def build_lookup_cache(self) -> None:
        """
        Precompute per-node attribute lists used by verify_path_concordance's hot path, so
        that function can do a single C-level list comparison instead of a per-node Python
        loop. Built lazily -- verify_path_concordance calls this itself on first use per
        anchor (checking self._node_ids is None) -- rather than eagerly for the whole
        anchor dictionary at load time: eager building scales with total anchor-dictionary
        size regardless of how many anchors any given run's reads actually touch, which is
        cheap on a small graph but a substantial, mostly-wasted cost on a large graph with a
        comparatively modest read set (confirmed: dominated main-process time on a
        981K-node graph with only 6,400 reads). These lists go stale if _nodes is mutated
        afterward -- anything that appends/inserts a node into an already-cached anchor
        must reset self._node_ids = None to force a rebuild.
        """
        self._node_ids = [n.id for n in self._nodes]
        self._node_ids_rev = self._node_ids[::-1]
        self._orientations = [n.orientation for n in self._nodes]
        _orientations_rev = self._orientations[::-1]
        # Precomputed negation so the discordant-orientation check is a direct list == too;
        # see verify_path_concordance for the derivation of why this is the right comparand.
        self._orientations_rev_negated = [not o for o in _orientations_rev]
        self._lengths = [n.length for n in self._nodes]
        self._lengths_rev = self._lengths[::-1]
        self._lengths_prefix = list(accumulate(self._lengths, initial=0))
        self._lengths_rev_prefix = list(accumulate(self._lengths_rev, initial=0))

    @profile
    def __len__(self):
        return len(self._nodes)
    
    @profile
    def __getitem__(self, position):
        return self._nodes[position]

    @profile
    def add(self, node):
        self._nodes.append(node)

    @profile
    def insert_node_through_extension(self, node, insert_left):
        if ((insert_left) and self.path_orientation) or ((not insert_left) and (not self.path_orientation)):
            self._nodes.insert(0, node)
        else:
            self._nodes.append(node)

    @profile
    def merge_anchor(self, new_anchor, insert_left=False) -> bool:
        if insert_left:
            if self._nodes[0].id != new_anchor[-1].id:
                return False
            self._nodes[:0] = new_anchor[:-1]
        else:
            if self._nodes[-1].id != new_anchor[0].id:
                return False
            self._nodes.extend(new_anchor[1:])
        return True

    @profile
    def flip_anchor(self):
        self._nodes = self._nodes[::-1]
        for node in self._nodes:
            node.orientation = not node.orientation

    @profile
    def add_snarl_id(self, snarl_id) -> None:
        self.snarl_id = snarl_id

    @profile
    def _snapshot(self) -> 'Anchor':
        """
        Lightweight pre-extension snapshot used by extend_and_merge_snarls, in place of
        copy.deepcopy(). During extension only two things are ever mutated in-place:
        _nodes (via insert_node_through_extension) and the inner lists of bp_matched_reads
        (via element assignment in _try_extension). Node objects and every other attribute
        are never mutated, so sharing them instead of deep-copying is safe and much cheaper.
        """
        snap = Anchor.__new__(Anchor)
        snap.__dict__.update(self.__dict__)
        snap._nodes = self._nodes[:]
        snap.bp_matched_reads = [r[:] for r in self.bp_matched_reads]
        return snap
    
    @profile
    def __repr__(self) -> str:
        anchor_str = ""
        for node in self._nodes:
            orientaiton = ">" if node.orientation else "<"
            anchor_str += orientaiton + str(node.id)
        return anchor_str

    @profile
    def bandage_representation(self) -> str:
        bandage_nodes_str = ""
        for node in self._nodes:
            bandage_nodes_str += "," + str(node.id)
        return bandage_nodes_str[1:]
    
    @profile
    def add_sequence(self) -> None:
        self.num_sequences += 1

    @profile
    def compute_snarl_boundary(self) -> None:
        """
        This function computes the start and end boundary of the anchor. 
        The start boundary is always the smaller node ID and end boundary is 
        always the larger node ID irrespective of the anchor orientation

        """
        if self._nodes[-1].id > self._nodes[0].id:
            self.snarl_start_node = self._nodes[0]
            self.snarl_end_node = self._nodes[-1]
        else:
            self.snarl_start_node = self._nodes[-1]
            self.snarl_end_node = self._nodes[0]


    @profile
    def compute_bp_length(self) -> int:
        """
        This function computes the size, in basepairs, of an anchor.
        First and last node (the ones at the boundaries) are just half the node, while the
        nodes in the middle are taken in full.

        Returns
        -------
        anchor_size: int
        The length in basepairs of the anchor
        """

        # self.baseparilength = (self._nodes[0].length // 2) + (self._nodes[-1].length // 2)
        # for node_handle in self._nodes[1:-1]:
        #     self.baseparilength += node_handle.length
        
        ## computing new basepairlength
        # self.basepairlength = (0 if self._nodes[0].length == 1 else 1) + (0 if self._nodes[-1].length == 1 else 1)
        
        self.basepairlength = self.bp_occupied_start_node + self.bp_occupied_end_node

        for node in self._nodes[1:-1]:
            self.basepairlength += node.length

        return 


    @profile
    def set_snarl_max_boundaries(self, start_node, start_node_bp_occupied, end_node, end_node_bp_occupied) -> None:
        self.snarl_max_left_boundary_node = start_node
        self.snarl_left_boundary_bp_occupied = start_node_bp_occupied
        self.snarl_max_right_boundary_node = end_node
        self.snarl_right_boundary_bp_occupied = end_node_bp_occupied


    @profile
    def get_sentinel_id(self) -> int:
        """
        This function takes a path traversing a snarl (candidate anchor) and returns the sentinel associated to it. The sentinel is not dependent on the anchor orientation.

        Returns
        -------
        sentinel: int
        the node_id of the sentinel node
        """
        # Determine the index based on the orientation of the first element
        if self._nodes[0].orientation:
            # If not reversed, count from the beginning
            index = (len(self._nodes) - 1) // 2
        else:
            # If reversed, count from the end
            index = -((len(self._nodes) + 1) // 2)

        return self._nodes[index].id
    

    @profile
    def get_sentinels(self) -> list:
        """
        This function returns the list of sentinel nodes for the candidate anchor.
        All nodes except the first and last node (the ones at the boundaries) are sentinel nodes.

        """
        return self._nodes[1:-1]
    

    @profile
    def compute_sentinel_bp_length(self) -> int:
        """
        This function computes the size, in basepairs, of sentinel nodes of an anchor.

        Returns
        -------
        sentinel_size: int
        The length in basepairs of all sentinel node(s)

        """

        self.sentinel_length = sum([sentinel_node.length for sentinel_node in self.get_sentinels()])
        
        return


    @profile
    def __eq__(self, other_anchor) -> bool:
        """
        This function verifies if two paths in the same snarls are equal, independently of the orientation.

        Parameters
        ----------
        path_1: list
        list of node handle objects of the first path
        path_2: list
        list of node handle objects of the first path

        Returns
        -------
        bool
        True if the paths are equal else False
        """

        # if different in length, don't need to bother
        if len(self._nodes) != len(other_anchor):
            return False

        pos_1 = pos_2 = 0
        orientation_concordance = self._nodes[pos_1].orientation == other_anchor[pos_2].orientation
        if not orientation_concordance:
            pos_1 = len(self._nodes) - 1

        for _ in range(len(self._nodes)):
            if self._nodes[pos_1].id != other_anchor[pos_2].id or orientation_concordance != (
                self._nodes[pos_1].orientation == other_anchor[pos_2].orientation
            ):
                return False
            pos_1 += 1 if orientation_concordance else -1
            pos_2 += 1
        return True
    
    @profile
    def is_preceding_anchor(self, other_anchor) -> bool:
        """
        This function checks if the current anchor is preceding the other anchor.
        """
        return min(self._nodes[0].id, self._nodes[-1].id) < min(other_anchor._nodes[0].id, other_anchor._nodes[-1].id)
    

    @profile
    def get_reference_paths(self):
        out_s = ""
        for el in self.reference_paths_covered:
            out_s += el + ','
        return out_s[:-1]


    @profile
    def get_bed(self):
        #CHROM CHROM_START CHROM_END NAME
        return f"{self.chromosome}\t{self.genomic_position}\t{self.genomic_position+self.baseparilength}\t{self.__repr__}"
