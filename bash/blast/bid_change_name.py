import sys

def parse_gff3(gff3_path):
    """Extracts mRNA to gene name mapping from a GFF3 file."""
    id_dict = {}
    with open(gff3_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes = parts[8]
            if feature_type == 'mRNA':
                try:
                    attr_dict = dict(item.split('=') for item in attributes.split(';') if '=' in item)
                    mRNA_id = attr_dict['ID'].split('.')[0]
                    gene_id = attr_dict['Parent']
                    id_dict[mRNA_id] = gene_id
                except (KeyError, IndexError, ValueError):
                    continue  # skip malformed attribute lines
    return id_dict

def convert_names(mapping_file, id_dict, output_file):
    with open(mapping_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                continue
            tokens = line.strip().split()
            mRNA_name = tokens[0].split('.')[0]
            if mRNA_name in id_dict:
                gene_name = id_dict[mRNA_name].replace('_', '-')
                ara_name = tokens[1].split('.')[0]
                new_line = [gene_name, ara_name] + tokens[2:]
                outfile.write('\t'.join(new_line) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python change_name.py <GFF3 file> <mapping file> <output file>")
        sys.exit(1)

    gff3_file = sys.argv[1]
    besthit_file = sys.argv[2]
    output_path = sys.argv[3]

    id_mapping = parse_gff3(gff3_file)
    convert_names(besthit_file, id_mapping, output_path)

