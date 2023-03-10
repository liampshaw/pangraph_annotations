# Aim:
# 
# INPUT
# pangraph input (json)
# file of annotations (for a strain that appears in pangraph)
# OUTPUT
# file of annotations with block and block positions added (?)

import argparse

def get_options():
    parser = argparse.ArgumentParser(description="Map existing annotations onto pancontigs in a graph")
    parser.add_argument("--pangraph", help="Pangraph input file (json)", type=str, required=True)
    parser.add_argument("--annotations", help="Annotation input file (gbk)", type=str, required=True)
    parser.add_argument("--output", help="Output file (?)", type=str, required=False)
    return parser.parse_args()

def extractPaths(paths):
    genome_dict = {}
    for path in paths:
        genome = path['name']
        positions = path['position'] # N.B. already 1-indexed
        blocks = [x['id'] for x in path['blocks']]
        strands = [strand(x['strand']) for x in path['blocks']]
        genome_dict[genome] = []
        for i, b in enumerate(blocks):
            if i==len(blocks)-1:
                offset = 0 # From documentation: If N is the number of nodes this list has N+1 entries. The last entry is the position of the right edge of the last block in the path.
            else:
                offset = -1
            genome_dict[genome].append([b, strands[i], positions[i], positions[i+1]+offset])
    return(genome_dict)



def main():
	args = get_options()

	with open(str(args.pangraph), 'r') as f:
        pangraph_json = json.load(f)

    genome_dict = extractPaths(pangraph_json['paths'])
    

if __name__=="__main__":
	main()