import os
import argparse
import re
from titan_to_gene import titan_to_gene

parser = argparse.ArgumentParser()
parser.add_argument("input_dir")
parser.add_argument("genes_file", help='path to the genes set file')
parser.add_argument("output_dir")
args = parser.parse_args()

for filename in os.listdir(args.input_dir):
    m = re.match(r"(.*)_cluster1.segs.txt", filename)
    if m:
        titan_to_gene('/'.join([args.input_dir, filename]), args.genes_file,
                      '/'.join((args.output_dir, m.group(1) + '_genes_cnv.csv')))
