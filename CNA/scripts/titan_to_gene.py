#
# Python script to map the segs.txt results from TitanCNA to gene set
#
# Usage:
# python ./titan_to_gene.py result_seg_file genes_set_file > output.txt
#

import pandas as pd
from collections import OrderedDict
import argparse
import progressbar


def titan_to_gene(result_file, genes_file, out=None, tel_genes=None):
    # load data
    data = pd.read_csv(result_file, delimiter='\t')
    genes = pd.read_csv(genes_file, delimiter='\t')

    # preprocessing
    interested_col = ['Sample', 'Chromosome', 'Start_Position.bp.', 'End_Position.bp.', 'Corrected_Call',
                      'Corrected_Copy_Number']
    data = data[interested_col]
    data['Chromosome'] = data['Chromosome']
    genes['Chrom'] = genes['Chrom'].map(lambda x: 'chr' + x)
    genes = genes[genes['Chrom'].isin(data['Chromosome'].unique())]
    genes = genes.rename(columns={'Chrom': 'Chromosome'})

    # merge genes and CNVs
    combine = pd.merge(genes, data, on=['Chromosome'])
    combine = combine.rename(columns={'Start': 'Gene_Start', 'End': 'Gene_End', 'Start_Position.bp.': 'CNV_Start',
                                      'End_Position.bp.': 'CNV_End'})
    # filter out unrelated rows
    combine = combine[(combine['Gene_End'] >= combine['CNV_Start']) & (combine['Gene_Start'] <= combine['CNV_End'])]
    keys = ['GeneID', 'GeneName', 'Chromosome', 'Gene_Start', 'Gene_End', 'Sample']

    # deal with multi-seg genes
    # full_cn: <seg1_copy_number>:<seg1_start>:<seg1_end>:<seg1_len>;<seg2_copy_number>:<seg2_start>:<seg2_end>2<seg2_len>
    res = list()
    count = 0
    bar = progressbar.ProgressBar(maxval=combine.shape[0],
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    for name, df in combine.groupby(keys):
        bar.update(count + 1)
        count += 1
        df = df.reset_index(drop=True)
        current = OrderedDict(zip(keys, name))
        current['max_cn'] = df['Corrected_Copy_Number'].max()  # report the max CN the gene has
        full_cn = list()
        for i in range(df.shape[0]):  # report each segment the gene covers
            full_cn.append(":".join([str(df['Corrected_Copy_Number'][i]), df['Corrected_Call'][i], str(df['CNV_Start'][i]),
                                     str(df['CNV_End'][i]), str(df['CNV_End'][i] - df['CNV_Start'][i])]))
        if current['max_cn'] > 2:
            current['gene_state'] = 1
        elif current['max_cn'] == 2:
            current['gene_state'] = 0
        else:
            current['gene_state'] = -1
        current['full_cn'] = ";".join(full_cn)

        res.append(current)
    bar.finish()
    res_df = pd.DataFrame(res).sort_values(['Chromosome', 'Gene_Start'])
    pd.options.display.max_colwidth = 10000

    # add telomere genes column
    if tel_genes:
        tel_genes = pd.read_csv(tel_genes)['Name']
        res_df['is_tel_related'] = res_df['GeneName'].isin(tel_genes).astype(int)

    # write to output
    if out:
        res_df.to_csv(out, index=False)
    else:
        print(res_df.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_file", help='path to the segment result file')
    parser.add_argument("genes_file", help='path to the genes set file')
    parser.add_argument("--out", help='output to a file', default=None)
    parser.add_argument("--tel-genes", action='store', dest='tel_genes',
                        help='the file which contains the names of telomere related genes', default=None)
    args = parser.parse_args()

    titan_to_gene(args.result_file, args.genes_file, args.out, args.tel_genes)
