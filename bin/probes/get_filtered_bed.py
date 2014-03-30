

import argparse
import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "--input",
            required=True,
            help="""The input BED file"""
        )
    parser.add_argument(
            "--output",
            required=True,
            help="""The output BED file"""
        )
    parser.add_argument(
            "--min-length",
            type=int,
            default=100,
            help="""The minimum length to filter""",
        )
    return parser.parse_args()

def main():
    args = get_args()
    with open(args.input, 'rU') as input:
        with open(args.output, 'w') as outf:
            for line in input:
                ls = line.strip().split('\t')
                if int(ls[2]) - int(ls[1]) >= args.min_length:
                    outf.write(line)

if __name__ == '__main__':
    main()
