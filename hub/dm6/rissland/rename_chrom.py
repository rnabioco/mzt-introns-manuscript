import sys
import argparse
import os

def parse_chromfile(fn):
    chroms = {}
    with open(fn) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            chroms[fields[0]] = fields[1]
    return chroms

def rename_chroms(fh, chrommap):
    
    for line in fh:
        fields = line.split("\t")

        if fields[0] in chrommap:
            fields[0] = chrommap[fields[0]]
            print("\t".join(fields), end = "")
        else:
            continue

def main():
    parser = argparse.ArgumentParser(description="""
    rename encode fastqs""")

    parser.add_argument('-bg',
                        '--bedgraph',
                        help ='bedgraph to rename chroms',
                        required = True)
    parser.add_argument('-c',
                        '--chrommap',
                        help ='renaming chromosome file ',
                        required = True)
    args=parser.parse_args()
    
    bg = args.bedgraph
    if bg == "-":
        fin = sys.stdin
    else:
        fin = open(bg, 'r')
    chrommap = args.chrommap

    chrs = parse_chromfile(chrommap)
    
    rename_chroms(fin, chrs)
    fin.close()

if __name__ == '__main__': main()
