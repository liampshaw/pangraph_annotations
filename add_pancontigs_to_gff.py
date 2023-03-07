import pandas as pd
import re
import argparse
from datetime import datetime

import pangraph_locator 
import pangraph_interface 

def get_options():
    parser = argparse.ArgumentParser(description="Add information on pancontig location to gff",
                                     prog="add_pancontigs_to_gff")
    parser.add_argument("--pangraph", 
        help="Input pangraph (JSON)", required=True)
    parser.add_argument("--input_gff", 
        help="Annotations (GFF)", required=True)
    parser.add_argument("--output_gff", 
        help="Output gff with pancontigs as attributes (GFF)", required=False, default="")
    parser.add_argument("--mode", choices=["keep_original", "make_new"],  
        help="Whether to keep original gff and add pancontig attributes (keep_original) or make a new gff wrt pancontigs (make_new)", required=False, default="keep_original")
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
        # Add pancontig info onto the gff
        self.new_gff = add_pancontigs_to_gff(self.pangraph_map, self.original_gff.gff)
        self.pancontig_gff = add_gff_to_pancontigs(self.pangraph_map, self.original_gff.gff)


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


def add_pancontig_info(pangraph_map, strain: str, gff_entry):
    """given a gff entry, adds panconting information"""
    blocks_for_gene_ids, blocks_for_gene_rel_pos, blocks_for_gene_occurrences = list(pangraph_map[strain].interval_to_blocks(gff_entry.start, gff_entry.end))
    new_gff_entries = []
    pancontigInfo = ",".join([blocks_for_gene_ids[i]+{True: "+", False: "-"}[blocks_for_gene_occurrences[i][2]]+
                        "_"+str(blocks_for_gene_occurrences[i][1]) for i in range(len(blocks_for_gene_ids))])
    # add a new attribute
    entry_attributes = gff_entry.attributes+";pancontigs="+pancontigInfo
    new_gff_entries.append(gffEntry([gff_entry.seqid, gff_entry.source, gff_entry.type, gff_entry.start, gff_entry.end, gff_entry.score, gff_entry.strand, gff_entry.phase, \
                                    entry_attributes]))
    return(new_gff_entries)


def add_pancontigs_to_gff(pangraph_map, original_gff):
    """adds pancontig information for each gff entry"""
    new_gff_list = []
    for gff_entry in original_gff:
        strain = gff_entry.seqid
        gff_partials = add_pancontig_info(pangraph_map, strain, gff_entry)
        for gff_partial in gff_partials:
            new_gff_list.append(gff_partial)
    return(GFF(new_gff_list))

def project_annotation_onto_pancontig(pangraph_map, strain: str, gff_entry):
    """given a gff entry, returns new gff entries mapped onto pancontigs (blocks) of pangraph"""
    blocks_for_gene_ids, blocks_for_gene_rel_pos, blocks_for_gene_occurrences = list(pangraph_map[strain].interval_to_blocks(gff_entry.start, gff_entry.end))
    new_gff_entries = []
    if len(blocks_for_gene_ids)==1: # if only one block, entry is not fragmented across multiple blocks
        entry_type = gff_entry.type # Inherit entry type
        block_strand = {True: '+', False: '-'}[blocks_for_gene_occurrences[0][2]] 
        block_occurrence = blocks_for_gene_occurrences[0][1]
        # NEED TO: convert start/end of gene to pancontig coordinates
        new_gff_entries.append(gffEntry([blocks_for_gene_ids[0], gff_entry.seqid, entry_type, blocks_for_gene_rel_pos[0][0], blocks_for_gene_rel_pos[0][1], gff_entry.score, gff_entry.strand, gff_entry.phase, \
                                    gff_entry.attributes+";pancontigID="+blocks_for_gene_ids[0]+";pancontigStrand="+block_strand+
                                    ";pancontigN="+str(block_occurrence)]))
    else: # otherwise, entry is fragmented across n>1 blocks
        parent_entry_id = re.sub("ID=", "", re.sub(";.*", "", gff_entry.attributes)) # parent entry ID from attributes
        parent_other_attributes = re.sub("^.*?;", "", gff_entry.attributes)
        entry_type = gff_entry.type
        gene_starting_block = blocks_for_gene_ids[0] # get the first and 
        #gene_starting_block_idx = 
        gene_ending_block = blocks_for_gene_ids[-1] # last block
        block_ids_strain = pangraph_map[strain].ids
        for i, b in enumerate(blocks_for_gene_ids): # we go through each block to output the gene fragments
            # We get the index of the block - accounts for duplicated blocks (a headache to understand but I think works ok)
            b_index = [i for i in range(len(block_ids_strain)) if block_ids_strain[i] in b][blocks_for_gene_occurrences[i][1]-1]
            block_strand = {True: '+', False: '-'}[blocks_for_gene_occurrences[i][2]] # get strand
            block_occurrence = blocks_for_gene_occurrences[i][1]
            b_start, b_end = pangraph_map[strain].b[b_index], pangraph_map[strain].e[b_index] # get block start and end
            #if b==gene_starting_block: # if starting block, return gene start
            if i==0:
                fragment_start, fragment_end = gff_entry.start, b_end
            if i==len(blocks_for_gene_ids)-1:
                fragment_start, fragment_end = b_start, gff_entry.end
            if i==0 or i==len(blocks_for_gene_ids)-1:
                fragment_type = entry_type
            else: # otherwise, return block start/end as the gene spans the block
                fragment_start, fragment_end = b_start, b_end
                fragment_type = entry_type
            if gff_entry.strand=="+":
                fragment_attributes = "ID="+parent_entry_id+"-fragment"+str(i+1)+";parent="+parent_entry_id+";"+parent_other_attributes+";pancontigID="+b+";pancontigStrand="+block_strand+";pancontigN="+str(block_occurrence)
                fragment_phase = calculate_phase(gff_entry.start, gff_entry.phase, fragment_start) #(((fragment_start-gff_entry.start) % 3) + gff_entry.phase) % 3 # correct the phase
            elif gff_entry.strand=="-": # reverse fragments if gene on negative strand
                fragment_attributes = "ID="+parent_entry_id+"-fragment"+str(len(blocks_for_gene_ids)-i)+";parent="+parent_entry_id+";"+parent_other_attributes+";pancontigID="+b+";pancontigStrand="+block_strand+";pancontigN="+str(block_occurrence)
                fragment_phase = calculate_phase(gff_entry.end, gff_entry.phase, fragment_end) # (((fragment_end-gff_entry.end) % 3) + gff_entry.phase) % 3 # correct the phase
            # Add position in block coordinates...
            #position_in_block_coordinates(fragment_start, b_start, b_end, )

            # Todo: - unclear what 'score' should be for a fragmented entry - inherited or not? Leaving blank ('.') for now
            new_gff_entries.append(gffEntry([blocks_for_gene_ids[i], gff_entry.seqid, \
                                    entry_type, \
                                    blocks_for_gene_rel_pos[i][0], blocks_for_gene_rel_pos[i][1], ".", gff_entry.strand, fragment_phase, fragment_attributes])) 
    return(new_gff_entries)



def add_gff_to_pancontigs(pangraph_map, original_gff):
    """adds a gff onto the pancontigs, making a note of fragmented genes in attributes"""
    new_gff_list = []
    for gff_entry in original_gff:
        strain = gff_entry.seqid
        gff_partials = project_annotation_onto_pancontig(pangraph_map, strain, gff_entry)
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
    additional_header_string = "#!pancontig information relative to "+str(args.pangraph)+" added on "+datetime.now().strftime("%m/%d/%Y, %H:%M:%S")+"\n"
    gff_header_string = gff_header(args.input_gff)+additional_header_string
    glued_gff = GraphGFF(args.pangraph, args.input_gff)
    if args.mode=="keep_original":
        output_gff_list = glued_gff.new_gff.gff_to_df().values.tolist()
    elif args.mode=="make_new":
        # IN PROGRESS
        # To do: add a proper header string with sequence regions as pancontigs
        # Need to output a new file of the pancontigs with *actual sequences* in the strain
        # so can inspect e.g. in IGV and match up
        output_gff_list = glued_gff.pancontig_gff.gff_to_df().values.tolist()
    if args.output_gff!="":
        write_gff(output_gff_list, args.output_gff, header_string = gff_header_string)
    else:
        print(gff_header_string)
        for entry in output_gff_list:
            print("\t".join([str(x) for x in entry]))





if __name__== "__main__":
    main()