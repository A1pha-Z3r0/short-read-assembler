# pylint: disable=trailing-whitespace

"""
This is the main entry point to test out a naive DBG based short read genome assembly.
"""

import argparse
from utils.helper import read_assembly


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Assemble short reads into contigs")
    parser.add_argument("input", help="Path to input .txt or .fasta file")
    parser.add_argument("output", help="Path to output .fasta file")
    parser.add_argument("-k", "--kmer", type=int, default=3, 
                        help="k-mer size for counting (default: 3)")

    return parser.parse_args()


def main():
    "This is the main function that has all the logic"
    args = parse_args()
    contigs = read_assembly(inp_file_path = args.input, 
                            out_file_path = args.output,
                            k = args.kmer)

  

    print(len(contigs))

if __name__ == "__main__":
    main()
