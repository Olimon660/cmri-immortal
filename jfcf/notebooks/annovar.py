import pandas as pd
from glob import glob
import re

df_list = []
#for fp in glob("/home/scai/processing_simoncai/tmm/tel_vs_alt/immortal_simon/annovar/raw/*"):
#    name = fp.replace("/home/scai/processing_simoncai/tmm/tel_vs_alt/immortal_simon/annovar/raw/","")
#    sample = re.search(r"(JF.*)_[INDEL|SNP].*", name).group(1)
#    df = pd.read_csv(fp, sep='\t')
#    df['sample'] = sample
#    df_list.append(df)
#df_all = pd.concat(df_list)
#df_all.to_csv("/home/scai/processing_simoncai/tmm/tel_vs_alt/immortal_simon/annovar/combined.csv", index=False)
df_all = pd.read_csv("/home/scai/processing_simoncai/tmm/tel_vs_alt/immortal_simon/annovar/combined.csv")
#df_all.to_pickle("/home/scai/processing_simoncai/tmm/tel_vs_alt/immortal_simon/annovar/combined.pkl")
df_all.to_hdf("./combined.hd", key='all')
