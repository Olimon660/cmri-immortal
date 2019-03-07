#
# Python script to map the segs.txt results from TitanCNA to gene set
#
# Usage:
# python ./genes_modelling.py /Users/scai/CMRI/immortal/results/all_genes \
#           --out /Users/scai/CMRI/immortal/results/ready_for_model/test.csv
#

import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("genes_cnv_dir", help='directory which contains the gene CNV results for each pair')
parser.add_argument("--out", help='output to a file')
args = parser.parse_args()

# Load in all gene results
all_genes = list()
labels = {'IIICF_402DE_D2': 0,
          'IIICF_d2': 0,
          'IIICF_E6_A1': 0,
          'IIICF_E6E7_A1': 0,
          'IIICF_E6E7_C4_post': 0,
          'IIICF_T_C3': 0,
          'LFS_05F_24_post': 0,
          'MeT_4A_post': 0,
          'IVG_BF_LXSN_post': 1,
          'GM847': 1,
          'IIICF_a2': 1,
          'IIICF_c': 1,
          'IIICF_E6_A2': 1,
          'IIICF_T_A6': 1,
          'IIICF_T_B1': 1,
          'IIICF_T_B3': 1,
          'IIICF_T_C4': 1,
          'JFCF_6_P_pLKO_5': 1,
          'JFCF_6_T_1_D': 1,
          'JFCF_6_T_1_L': 1,
          'JFCF_6_T_1_M': 1,
          'JFCF_6_T_1_P_ALT': 1,
          'JFCF_6_T_1_R': 1,
          'JFCF_6_T_1J_1_3C': 1,
          'JFCF_6_T_1J_11E': 1,
          'JFCF_6_T_5K': 1,
          'JFCF_6_T_1J_6B': 1,
          'VA13': 1,
          'JFCF_6_T_1_C': 2,
          'JFCF_6_T_1_F': 2,
          'JFCF_6_T_1_G': 2,
          'JFCF_6_T_1_H': 2,
          'JFCF_6_T_1_P_TEL': 2,
          'JFCF_6_T_1J_11C': 2,
          'JFCF_6_T_2H': 2,
          'JFCF_6_T_1_Q': 2
          }
for file in os.listdir(args.genes_cnv_dir):
    df = pd.read_csv(args.genes_cnv_dir + '/' + file)
    sample = df['Sample'][0]
    df = df[['GeneName', 'gene_state']].rename(columns={'gene_state': sample})
    all_genes.append(df)

# Merge results from all samples
all_genes_df = all_genes[0]
for i in range(1, len(all_genes)):
    all_genes_df = pd.merge(all_genes_df, all_genes[i], how='outer', on='GeneName', sort=False)
all_genes_df = all_genes_df.fillna(0)

# Wrangle the data so the format is ready for sklearn
all_genes_df_t = all_genes_df.transpose()
all_genes_df_t.columns = all_genes_df_t.iloc[0, :]
all_genes_df_t = all_genes_df_t.iloc[1:, :]
all_genes_df_t = all_genes_df_t.reset_index()
all_genes_df_t = all_genes_df_t.rename(columns={'index': 'sample'})
print(all_genes_df_t.shape)
all_genes_df_t['y'] = all_genes_df_t['sample'].map(lambda x:labels[x])
all_genes_df_t = all_genes_df_t.sort_values(['y'])

# Output
if args.out:
    all_genes_df_t.to_csv(args.out, index=False)
else:
    print(all_genes_df_t.to_string(index=False))
