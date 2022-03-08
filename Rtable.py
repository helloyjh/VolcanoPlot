import os
import statsmodels
import statsmodels.api as sm
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Counter')
parser.add_argument('--rank', type=str, default='O')
parser.add_argument('--m', type=str, default='DESeq2')
parser.add_argument('--v', type=str, default='1.05')
parser.add_argument('--p', type=int, default=20)
args = parser.parse_args()

if args.m == 'DESeq2':
    xy = ['log2FoldChange', 'padj']
elif args.m == 'edgeR':
   xy = ['logFC', 'FDR']

p_path = 'R_result/{}/p_rank_{}.tsv'.format(args.m, args.rank)

data = pd.read_csv(p_path, sep='\t')

data_adjp = data.sort_values(by=xy[1])
data_FC = data.sort_values(by=xy[0])
data_adjp = data_adjp[:args.p][['Unnamed: 0', xy[0], xy[1]]]
data_FCL = data_FC[:int(args.p/2)][['Unnamed: 0', xy[0], xy[1]]]
data_FCR = data_FC[-int(args.p/2):][['Unnamed: 0', xy[0], xy[1]]]
data_FC = pd.concat([data_FCL, data_FCR])

os.makedirs('R_result/{}'.format(args.m), exist_ok=True)
data_FC.to_csv('R_result/{}/FC_p{}_rank_{}.tsv'.format(args.m, args.p, args.rank), index=False, sep='\t')
data_adjp.to_csv('R_result/{}/adjp_p{}_rank_{}.tsv'.format(args.m, args.p, args.rank), index=False, sep='\t')
