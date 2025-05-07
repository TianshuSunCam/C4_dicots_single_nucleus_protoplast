import sys

def get_best_hit(input_file, output_file):
    num_line = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if not line:  # skip empty lines
                continue
            if line.startswith('#'):
                num_line = 0
            else:
                num_line += 1
                if num_line == 1:
                    outfile.write(line + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_besthit.py <input_file> <output_file>")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
    get_best_hit(input_filename, output_filename)
