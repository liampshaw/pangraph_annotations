# Author: Marco Molari
# Source: https://github.com/mmolari/pangenome_evolution/blob/main/src/py_utils/pangraph_interface.py

# Wrapper to import in python the results of the Pangraph pipeline.

import numpy as np
import json
import pandas as pd

from collections import Counter
import pangraph_alignment as pga


def run_pangraph(align, output, compressed=False):
    """TODO: implement a wrapper to run pangraph from python? Possibly save
    the results in compressed format?
    """
    pass


class Pangraph:
    """Wrapper class to load and interact with the output of the Pangraph pipeline.
    The class has two main attributes:
    - `paths` : each strain has
    - `blocks` :
    """

    def __init__(self, pan_json):
        """Python calss to load the output of the Pangraph pipeline.

        Args:
            pan_json (dict): content of the .json file produced by pangraph.
        """
        self.paths = PathCollection(pan_json["paths"])
        self.blocks = BlockCollection(pan_json["blocks"])

    @staticmethod
    def load_json(filename):
        """Creates a Pangraph object by loading it from the .json file.

        Args:
            load_json (str): .json file to be loaded.

        Returns:
            Pangraph: the Pangraph object containing the results of the pipeline.
        """

        isjson = str(filename).endswith(".json")
        if not isjson:
            raise Exception(f"the input file {filename} should be in .json format")

        with open(filename, "r") as f:
            pan_json = json.load(f)
        pan = Pangraph(pan_json)
        return pan

    def strains(self):
        """Return lists of strain names"""
        return self.paths.ids_copy()

    def block_ids(self):
        """Returns the list of block ids"""
        return self.blocks.ids_copy()

    def to_json(self):
         # TODO: implement saving to json format after manipulation.
         # Maintain the same input format
         pass

    def to_paths_dict(self):
        """Generates a compressed representation of paths as simply lists of
        block ids. Returns a dictionary strain name -> list of block ids.
        """
        return self.paths.to_block_dict()

    def to_blockcount_df(self):
        """Returns a dataframe whose rows are strain names, and columns are block
        names. Values indicate the number of times a block is present. This can
        also be used to build a presence / absence matrix."""
        block_counters = {}
        for path in self.paths:
            block_counters[path.name] = Counter(path.block_ids)
        return pd.DataFrame(block_counters).fillna(0).T

    def to_blockstats_df(self):
        """Returns a dataframe containing statistics about blocks distribution.
        The index of the dataframe are block ids, and the columns are:
        - count: n. times the block occurs
        - n. strains: number of strains in which the block is observed
        - duplicated: whether the block is duplicated in at least one strain
        - len: average block length from pangraph.
        - core: whether a gene occurrs exactly once per strain
        """
        block_counter = Counter()
        str_counter = Counter()

        for path in self.paths:
            block_counter.update(path.block_ids)
            str_counter.update(np.unique(path.block_ids))

        lengths = {block.id: len(block.sequence) for block in self.blocks}

        df = pd.DataFrame(
            {"count": block_counter, "n. strains": str_counter, "len": lengths}
        )
        df["duplicated"] = df["count"] > df["n. strains"]
        df["core"] = (df["n. strains"] == len(self.paths)) & (df["duplicated"] == False)
        return df


class IndexedCollection:
    """This class is used to implement smart indexing of a list of blocks or paths.
    The class has three elements:
    - ids: ordered list of item ids (string)
    - list: ordered list of items
    - id_to_pos: dictionary mapping ids to their position in the list.

    An element of the class can be indexed in these different ways:
    - through its id (string)
    - through its position on the list
    - though a list of ids
    - through a list of positions
    - though a boolean mask on the list

    This is handled by the __getitem__ function.

    Moreover the object can be cast to iterators, in which case an iterator over
    the list items is returned.

    The object's `len` is the lentgth of the list.
    """

    def __init__(self, ids, items):
        self.ids = np.array(ids)
        self.list = np.array(items)
        self.id_to_pos = {id: n for n, id in enumerate(ids)}

    def __iter__(self):
        return iter(self.list)

    def __len__(self):
        return len(self.list)

    def __getitem__(self, idx):

        # if indexed by block id
        if isinstance(idx, str):
            pos = self.id_to_pos[idx]
            return self.list[pos]

        # if indexed by integer
        if isinstance(idx, (int, np.integer)):
            return self.list[idx]

        # if indexed by list or numpy array
        if isinstance(idx, (list, np.ndarray)):
            # if list is empty return empty list
            if len(idx) == 0:
                return []

            idx0 = idx[0]
            # if the type is integer, return corresponding items
            if isinstance(idx0, (int, np.integer)) and not isinstance(idx0, bool):
                return self.list[idx]

            # if the type is string, return corresponding ids
            if isinstance(idx0, str):
                return [self.__getitem__(idx_i) for idx_i in idx]

            # if the type is bool (a mask)
            if isinstance(idx0, np.bool_):
                return self.list[idx]

        # if no condition is matched, then raise an error
        message = """
        the index object passed does not match any of the allowed types:
        - integer or string
        - list of integers or strings
        - boolean numpy array (mask)
        """
        raise TypeError(message)

    def ids_copy(self):
        return self.ids.copy()


class BlockCollection(IndexedCollection):
    """Collection of all blocks. Inherits from IndexedCollection to allow for
    smart indexing of blocks.
    """

    def __init__(self, pan_blocks):
        ids = [block["id"] for block in pan_blocks]
        items = [Block(block) for block in pan_blocks]
        IndexedCollection.__init__(self, ids, items)


class PathCollection(IndexedCollection):
    """Collection of all paths. Inherits from IndexedCollection to allow for
    smart indexing of paths.
    """

    def __init__(self, pan_paths):
        ids = [path["name"] for path in pan_paths]
        items = [Path(path) for path in pan_paths]
        IndexedCollection.__init__(self, ids, items)

    def to_block_dict(self):
        return {path.name: path.block_ids.copy() for path in self}


class Block:
    """Python wrapper for pangraph block object. It has attributes:
    - id (str): block id
    - sequence (str): sequence of the block
    - gaps (list): list of gap occurrences in the alignment, as returned by
        pangraph
    - alignment (pan_alignment): object containing the information provided by
        pangraph that can be used to build alignments. See `pan_alignment`
        class for details.
    """

    def __init__(self, pan_block):
        self.id = pan_block["id"]
        self.sequence = pan_block["sequence"]
        self.alignment = pga.pan_alignment(pan_block)

    def __len__(self):
        """Length of the sequence in base-pairs."""
        return len(self.sequence)

    def __str__(self):
        return f"block {self.id}, consensus len {len(self.sequence)/1000} kbp, {self.depth()} occurrences."

    def depth(self):
        """How many occurrences of the block are present"""
        return len(self.alignment.occs)

    def frequency(self):
        """In how many unique strains is the block present"""
        return len(self.strains())

    #def strains(self):
    #    """Which (unique) strains have the block."""
    #    return np.unique([w[0] for w in self.alignment.occs])

    def is_duplicated(self):
        """Check whether a block is duplicated in any strain"""
        return self.depth() > self.frequency()


class Path:
    """Python wrapper for pangraph block object. It has attributes:
    - name (str): strain name
    - offset (int): offset of the first block in path?
    - circular (bool): whether the path object is circular

    then there are three attributes that are lists. Each index represent a
    particular block occurrence. These are
    - block_ids (str) : list of block ids
    - block_nums (int) : list of block occurrence (blocks can occurr more than
        once per strain)
    - block_strains (str) : list of strains in which blocks occurr
    - block_strands (bool) : whether the block occurrs on direct or reverse strand.
    """

    def __init__(self, pan_path):
        self.name = pan_path["name"]
        self.offset = pan_path["offset"]
        self.circular = pan_path["circular"]
        blocks = pan_path["blocks"]
        self.block_ids = np.array([block["id"] for block in blocks])
        self.block_nums = np.array([block["number"] for block in blocks])
        self.block_strains = np.array([block["name"] for block in blocks])
        self.block_strands = np.array([block["strand"] for block in blocks])

    def __len__(self):
        return len(self.block_ids)

    def __str__(self):
        return f"path {self.name}, n. blocks = {len(self.block_ids)}"