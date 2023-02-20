import pandas as pd
import re
import argparse
from datetime import datetime

import pangraph_locator 
import pangraph_interface 

def get_options():
    parser = argparse.ArgumentParser(description="Add gff annotations onto a pangraph object",
                                     prog="glue_gff")
    parser.add_argument("--pangraph", help="Input pangraph (JSON)", required=True)
    parser.add_argument("--input_gff", help="Annotations (GFF)", required=True)
    parser.add_argument("--output_gff", help="Output gff with pancontigs as attributes (GFF)", required=False, default="")
    return parser.parse_args()

class gffEntry:
    """Class for a single gff entry"""
    def __init__(self, gff_list):
        """Takes a list of 9 objects and creates the gff entry"""
        self.seqid = gff_list[0]
        self.source = gff_list[1]
        self.type = gff_list[2]
        self.start = int(gff_list[3])
        self.end = int(gff_list[4])
        self.score = gff_list[5]
        self.strand = gff_list[6]
        self.phase = gff_list[7] # phase is '.' for everything else (e.g. 'region') so should not int
        self.attributes = gff_list[8]

class GFF:
    """A data holder for gff."""
    def __init__(self, gff_file_or_list):
        if type(gff_file_or_list) is list:
            self.gff = gff_file_or_list
        elif type(gff_file_or_list) is str:
            self.gff = load_gff(gff_file_or_list)    
    def gff_to_df(self):
        """returns a gff as pandas dataframe"""
        gff_df = pd.DataFrame([vars(x) for x in self.gff])
        return gff_df

class GraphGFF:
    """A data holder for original gff for a strain and the new gff mapped onto a pangraph."""
    def __init__(self, pangraph_file, gff_file):
        self.original_gff = GFF(gff_file)
        self.pangraph = pangraph_interface.Pangraph.load_json(pangraph_file)
        # Locator
        loc = pangraph_locator.Locator(self.pangraph)
        # Get map of the pangraph
        self.pangraph_map = pangraph_locator.build_map(self.pangraph.paths, self.pangraph.blocks)
        # 'Glue' the gff to the graph
        self.new_gff = glue_gff(self.pangraph_map, self.original_gff.gff)


def file_is_gff(gff_file):
    """Checks if a file has the expected gff header"""
    with open(gff_file, 'r') as f:
        try:
            gff_comment = f.readline().split(' ')[0]
        except:
            return(False)
        if gff_comment!="##gff-version":
            return(False)
        else:
            return(True)

def load_gff(gff_file):
    """Loads in a gff file"""
    if file_is_gff(gff_file):
        gff_entry_list = []
        with open(gff_file, "r") as f:
            for line in f.readlines():
                if line.startswith("#"):
                    pass
                if line.startswith("##FASTA"): # Don't read in fasta components
                    break
                else:
                    line = line.strip("\n").split("\t")
                    if len(line)==9:
                        gff_entry_list.append(line)
    gff = [gffEntry(x) for x in gff_entry_list]
    return(gff)

def gff_header(gff_file):
    """Extracts the header of a gff_file"""
    if file_is_gff(gff_file):
        gff_header_string = ""
        with open(gff_file, "r") as f:
            for line in f.readlines():
                if line.startswith("#"):
                    gff_header_string += line
                else:
                    break
    return(gff_header_string)

def rename_id(gff):
    """renames features using the seqid, which is more interpretable later on 
    ID=1_1 -> ID=seqid_1
    """
    new_gff = list(gff)
    for i in range(0, len(gff)):
        seqid = str(gff[i].seqid)
        new_attribute = re.sub("ID=.*_(\d+);", "ID="+seqid+"_\\1;", gff[i].attributes)
        new_gff[i].attributes = str(new_attribute)
    return(new_gff)


def create_new_gff_entries(pangraph_map, strain: str, gff_entry):
    """given a gff entry, returns new gff entries mapped onto blocks of pangraph"""
    blocks_for_gene_ids, blocks_for_gene_rel_pos, blocks_for_gene_occurrences = list(pangraph_map[strain].interval_to_blocks(gff_entry.start, gff_entry.end))
    new_gff_entries = []
    pancontigInfo = ",".join([blocks_for_gene_ids[i]+{True: "+", False: "-"}[blocks_for_gene_occurrences[i][2]]+
                        "_"+str(blocks_for_gene_occurrences[i][1]) for i in range(len(blocks_for_gene_ids))])
    # add a new attribute
    entry_attributes = gff_entry.attributes+";pancontigs="+pancontigInfo
    new_gff_entries.append(gffEntry([gff_entry.seqid, gff_entry.source, gff_entry.type, gff_entry.start, gff_entry.end, gff_entry.score, gff_entry.strand, gff_entry.phase, \
                                    entry_attributes]))
    return(new_gff_entries)


def glue_gff(pangraph_map, original_gff):
    """adds a gff onto the blocks of the graph, fragmenting entries if necessary across blocks"""
    new_gff_list = []
    for gff_entry in original_gff:
        strain = gff_entry.seqid
        gff_partials = create_new_gff_entries(pangraph_map, strain, gff_entry)
        for gff_partial in gff_partials:
            new_gff_list.append(gff_partial)
    return(GFF(new_gff_list))

def calculate_phase(start, initial_phase, position):
    if initial_phase==".":
        return(".")
    else:
        relative_phase_to_start = abs(position-start) % 3
        absolute_phase = (relative_phase_to_start + int(initial_phase)) % 3
        return(absolute_phase)

def write_gff(gff_list, gff_file, header_string="##gff-version 3\n"): 
    """writes a gff to file"""
    with open(gff_file, "w") as f:
            f.write(header_string)
            for entry in gff_list:
                f.write("\t".join([str(x) for x in entry])+"\n")

def main():
    args = get_options()
    additional_header_string = "#!pancontig information relative to "+str(args.pangraph)+" added by glue_gff on "+datetime.now().strftime("%m/%d/%Y, %H:%M:%S")+"\n"
    gff_header_string = gff_header(args.input_gff)+additional_header_string
    glued_gff = GraphGFF(args.pangraph, args.input_gff)
    output_gff_list = glued_gff.new_gff.gff_to_df().values.tolist()
    if args.output_gff!="":
        write_gff(output_gff_list, args.output_gff, header_string = gff_header_string)
    else:
        print(gff_header_string)
        for entry in output_gff_list:
            print("\t".join([str(x) for x in entry]))



if __name__== "__main__":
    main()