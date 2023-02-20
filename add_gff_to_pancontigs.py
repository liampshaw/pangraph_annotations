
# We go the other way
# And take as input:
# - fasta file of pancontigs
# - gff file of annotations onto a genome used to build those pancontigs
# Output:
# - gff file for the annotations projected onto the pancontigs 
# I suppose it has to be with respect to a single strain? 
# Or we can use 'source' information for the strain...

##sequence-region NZ_CP055986.1 1 206301

def get_fasta_info(fasta_file):
    """Reads a fasta file and returns dict of id and length of sequence"""
    seqs = {}
    current_seq_id = ""
    for line in open(fasta_file, "r").readlines():
        if line.startswith(">"):
            seq_id = line.strip("\n")[1:].split(" ")[0]
            seqs[seq_id] = 0
            current_seq_id = seq_id
        elif current_seq_id != "":
                seqs[current_seq_id] += len(line.strip("\n"))
    return(seqs)


def gff_regions(fasta_file):
    """Turns a fasta file into a set of regions for gff output"""
    



#def annotations_on_blocks(gff_list, pangraph):
# GENOME,BLOCK_ID,BLOCK_OCCURRENCE,BLOCK_STRAND,ANNOTATION_TYPE,
def annotations_on_blocks(pangraph_map, strain: str, gff_entry):
    """given a gff entry, annotations on blocks in a pseudo-gff format
    which doesn't meet gff requirements but may be useful"""
    blocks_for_gene_ids, blocks_for_gene_rel_pos, blocks_for_gene_occurrences = list(pangraph_map[strain].interval_to_blocks(gff_entry.start, gff_entry.end))
    block_annotations = []
    if len(blocks_for_gene_ids)==1: # if only one block, entry is not fragmented across multiple blocks
        entry_type = gff_entry.type # Inherit entry type
        block_strand = {True: '+', False: '-'}[blocks_for_gene_occurrences[0][2]] 
        block_occurrence = blocks_for_gene_occurrences[0][1]
        block_annotations.append(gffEntry([gff_entry.seqid, gff_entry.source, entry_type, gff_entry.start, gff_entry.end, gff_entry.score, gff_entry.strand, gff_entry.phase, \
                                    gff_entry.attributes+";pancontigID="+blocks_for_gene_ids[0]+";pancontigStrand="+block_strand+
                                    ";pancontigN="+str(block_occurrence)]))
    else: # otherwise, entry is fragmented across n>1 blocks
        parent_entry_id = re.sub("ID=", "", re.sub(";.*", "", gff_entry.attributes)) # parent entry ID from attributes
        parent_other_attributes = re.sub("^.*?;", "", gff_entry.attributes)
        entry_type = gff_entry.type+"-fragment" # Entry type is partial
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
                fragment_type = entry_type+"-edge"
            else: # otherwise, return block start/end as the gene spans the block
                fragment_start, fragment_end = b_start, b_end
                fragment_type = entry_type+"-middle"
            if gff_entry.strand=="+":
                fragment_attributes = "ID="+parent_entry_id+"-fragment"+str(i+1)+";parent="+parent_entry_id+";"+parent_other_attributes+";pancontigID="+b+";pancontigStrand="+block_strand+";pancontigN="+str(block_occurrence)
                fragment_phase = calculate_phase(gff_entry.start, gff_entry.phase, fragment_start) #(((fragment_start-gff_entry.start) % 3) + gff_entry.phase) % 3 # correct the phase
            elif gff_entry.strand=="-": # reverse fragments if gene on negative strand
                fragment_attributes = "ID="+parent_entry_id+"-fragment"+str(len(blocks_for_gene_ids)-i)+";parent="+parent_entry_id+";"+parent_other_attributes+";pancontigID="+b+";pancontigStrand="+block_strand+";pancontigN="+str(block_occurrence)
                fragment_phase = calculate_phase(gff_entry.end, gff_entry.phase, fragment_end) # (((fragment_end-gff_entry.end) % 3) + gff_entry.phase) % 3 # correct the phase
            # Add position in block coordinates...
            #position_in_block_coordinates(fragment_start, b_start, b_end, )

            # Todo: - unclear what 'score' should be for a fragmented entry - inherited or not? Leaving blank ('.') for now
            block_annotations.append(gffEntry([gff_entry.seqid, gff_entry.source, fragment_type, fragment_start, fragment_end, ".", gff_entry.strand, fragment_phase, fragment_attributes])) 
    return(block_annotations)
