import pandas as pd
import numpy as np
import os
import argparse
parser = argparse.ArgumentParser()

# Inputs --path:path to the folder of input PSM file ; --i:input file name ; --o: output file name ; --type: search algorithum used, for global proteomics use "Mascot", for MHC use "Sequest"

parser.add_argument('--path', type=str, required=True)
parser.add_argument('--i', type=str, required=True)
parser.add_argument('--o', type=str,required=True)
parser.add_argument('--type', type=str,required=True)
parser.add_argument('--q', type=float,required=True)

args = parser.parse_args()
os.chdir(args.path)
psm_table=pd.read_csv(args.i, sep='\t',engine='python')

def PSM_processing_Mascot(psm_table, Ion_score, search_rank):
    peptide_length=np.zeros(len(psm_table['Sequence']))
    for i in range(0,len(psm_table['Sequence'])):
        peptide_length[i]=len(list(psm_table['Sequence'][i]))
    psm_table['Seq_len']=peptide_length
    # Filter out high quality psms
    high_rank=psm_table[psm_table['Search Engine Rank']==search_rank]
    high_qvalue=high_rank[high_rank['Ions Score']>=Ion_score]
    return(high_qvalue)
def PSM_processing(psm_table,min_length, max_length, qvalue):
    peptide_length=np.zeros(len(psm_table['Sequence']))
    for i in range(0,len(psm_table['Sequence'])):
        peptide_length[i]=len(list(psm_table['Sequence'][i]))
    psm_table['Seq_len']=peptide_length
    # Filter out high quality psms
    high_qvalue=psm_table[psm_table['Percolator q-Value']<=qvalue]
    long_length=high_qvalue[high_qvalue['Seq_len']>=min_length]
    good_length=long_length[long_length['Seq_len']<=max_length]
    return(good_length)
#Default parameter 
min_length=8
max_length=12
max_length2=30
qvalue=args.q
search_rank=1
if (args.type=="Mascot"):
    good_length=PSM_processing_Mascot(psm_table,20, 1)
else:
    good_length=PSM_processing(psm_table,min_length, max_length, qvalue)
good_length.to_csv(args.o)
print('Done')