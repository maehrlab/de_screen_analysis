from __future__ import print_function
import numpy as np
import pandas as pd
import sys
import os
import re
PATH_TO_MIMOSCA=os.path.expanduser("~/Desktop/programs/MIMOSCA-master")
sys.path.append(PATH_TO_MIMOSCA)
from mimosca import read_10x, run_model_bycol, Zgenes, tp10k_transform, make_shufs, fdr_colwise_coefs, DE2GO

# =============== Run MIMOSCA on our DE perturb-seq data ============== #
# This script is a simple, hardwired way to run MIMOSCA on our samples.
# It loads the data from 'data/<dirname>', where dirname is usually 'DE_TERA_c0'.
# It expects 10X formatted data, plus a tsv file with metadata and the same cell names.
# It performs within-cell normalization, log1p, and within-gene standardization.
# It pulls the column 'highest_expressed_guide' from the metadata.
# It runs MIMOSCA with default params and EM-like correction of guide assignments.
# Results go into 'results/<dirname>', including coeffs and corrected assignments but not residuals.
# ===================================================================== #


def loadData( dirname, test ):
    print("Reading data.")
    Y = read_10x(os.path.join('data', dirname))
    Y.index = [re.sub('PLACEHOLDER_FOR_ENSG_', '', gene) for gene in Y.index]
    X = pd.read_table(os.path.join('data', dirname, 'metadata.tsv'))
    if test:
        X = X.head(100)
        Y = Y.transpose().head(100).transpose()
    assert all(X.index==Y.columns)
    return X, Y

def modelMatrix( my_factor ):
    print("Forming model matrix.")
    assert isinstance(my_factor, pd.core.series.Series)
    assert isinstance(my_factor.iloc[0], str)
    rownames = my_factor.index.copy()
    levels = list(set(my_factor))
    mod_mat = my_factor[:,None] == levels
    mod_mat = mod_mat.astype(int)
    mod_mat = pd.DataFrame(mod_mat)
    mod_mat.columns = levels
    mod_mat.index = rownames.copy()
    return mod_mat

def make_rp(dirname, test):
    if test:
        results_path_test_indicator = "_test=" + str(test)
    else:
        results_path_test_indicator = ""
    rp = dirname + results_path_test_indicator
    rp = os.path.join("results", rp)
    try:
        os.makedirs(rp)
    except:
        pass
    return rp

def load_run_save( dirname, test ):

    # Data loading
    X, Y = loadData( dirname, test )
    print(Y.iloc[0:3, 0:3])
    print(X.iloc[0:3, 0:3])
    Y = np.log2( 1 + tp10k_transform( Y ) )
    high_expressed_genes = list(np.sum(Y, axis = 1).sort_values(ascending = False)[0:1000].index)
    Y = Zgenes(Y)
    # Study interactions between guides and replicates: coeffs should be consistent.
    guide_and_rep = pd.Series([g + "_" + r for g,r in zip(X['highest_expressed_guide'], X['orig.ident'])])
    guide_and_rep.index = X.index.copy()
    mod_mat = modelMatrix(guide_and_rep)
    trt_guides = [guide_name for guide_name in mod_mat.columns if not re.match("Scramble", guide_name)]
    mod_mat = mod_mat[trt_guides]
    assert all(mod_mat.index == Y.columns)
    # MIMOSCA estimation
    print("Running MIMOSCA.")
    Be,X_adjust,RES_out = run_model_bycol( Y.transpose(),
                                           mod_mat,
                                           EM_cols=list(mod_mat.columns),
                                           modalpha=0.005,verbose=0 )
    print("Saving output.")
    rp = make_rp(dirname, test)
    Be.to_csv(      path_or_buf=os.path.join(rp, "coeffs.csv"))
    X_adjust.to_csv(path_or_buf=os.path.join(rp, "indicators_adjusted.csv"))

    # Permutation testing
    print("Running MIMOSCA permutation testing.")
    Be_shuffs = make_shufs(X = mod_mat, Xother = None, Y = Y.transpose(), shufnum=3)
    pvals = fdr_colwise_coefs(Be, Be_shuffs)
    pvals.to_csv(  path_or_buf = os.path.join(rp, "pvals.csv"))

    # GO terms
    print("Running MIMOSCA GO term analysis.")
    if test:
        pvals = pvals.head(200).copy()
        pvals = pvals.transpose().head(200).transpose().copy()
    df_bigGO = DE2GO(df_p=pvals, background=pvals.index, sig_thresh=np.log10(0.2), fdr_thresh=0.2, species='human')
    df_bigGO.to_csv(path_or_buf=os.path.join(rp, "GO.csv"))

    return

load_run_save('DE_TERA', test = True)
load_run_save('DE_TERA', test = False)
load_run_save('DE_TERA_c0', test = True)
load_run_save('DE_TERA_c0', test = False)

