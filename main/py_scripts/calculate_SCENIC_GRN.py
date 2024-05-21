#modifying functions
import pandas as pd
import pickle
import ast
import numpy as np 
import os

np.object = object

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

#from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from dask.diagnostics import ProgressBar
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

#define function to calculate regulons from consensus GRN ADJ and CSV files
#note that the dbs and motif_annotations are the same for all contrasts
def calculate_regulons(C_ADJ_FNAME, CSV_FNAME, dbs, motif_annotations, auc_threshold=0.05, nes_threshold=3.0):
    print('Calculating regulons')
    print('With parameters: auc_threshold =', auc_threshold, 'nes_threshold =', nes_threshold)
    adj_GRN = pd.read_csv(C_ADJ_FNAME, sep='\t')
    csv_GEX = pd.read_csv(CSV_FNAME, index_col=0)
    
    modules = list(modules_from_adjacencies(adj_GRN, csv_GEX))
    df = prune2df([dbs], modules, motif_annotations, auc_threshold=auc_threshold, nes_threshold=nes_threshold)
    regulons = df2regulons(df)
    
    return regulons, modules

#calculate SCENIC GRN from regulons object and save to results folder
def calculate_SCENIC_GRN(regulons): 
    SCENIC_GRN = pd.DataFrame(columns=['TF', 'target', 'importance'])
    for regulon in regulons:
        regulon_name = regulon.name
        if 'Regulon for' in regulon_name:
            regulon_name = regulon_name.split('Regulon for ')[1]
        gene_importance_pairs = pd.DataFrame((regulon.gene2weight.items()))
        gene_importance_pairs['TF'] = regulon_name
        gene_importance_pairs.columns = ['target', 'importance', 'TF']
        gene_importance_pairs = gene_importance_pairs[['TF', 'target', 'importance']]
        SCENIC_GRN = pd.concat([SCENIC_GRN, gene_importance_pairs])
    return SCENIC_GRN

#SCENIC_GRN = calculate_SCENIC_GRN(regulons)
#SCENIC_GRN


db_fpath = "/sc/arion/projects/roussp01a/collin/scenic-kaiyi/data/db_f"
f_db = os.path.join(db_fpath, 'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')
motif_annotations = os.path.join(db_fpath, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')
tf_list = os.path.join(db_fpath, 'allTFs_hg38.txt')

db_fnames = f_db
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = RankingDatabase(f_db, 'hg38_10kbup_10kbdown_fulltx_v10_clust')

