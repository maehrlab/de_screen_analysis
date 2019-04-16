from __future__ import print_function
import pandas as pd
import os
import math
from goatools.associations import read_ncbi_gene2go
import re

def invert_dol_nonunique(d):
    newdict = {}
    for k in d:
        for v in d[k]:
            newdict.setdefault(v, []).append(k)
    return newdict

PATH_TO_MIMOSCA = os.path.expanduser(os.path.join("~", "Desktop", "programs", "MIMOSCA-master"))
path2db   = os.path.expanduser(os.path.join(PATH_TO_MIMOSCA, "common_files", "data")) + os.path.sep
path2Xref = os.path.expanduser(os.path.join(PATH_TO_MIMOSCA, "common_files", "data")) + os.path.sep

geneid2gos = read_ncbi_gene2go(path2db + "gene2go", taxids=[9606])
go2idx = invert_dol_nonunique(geneid2gos)

Xtable = pd.read_csv(os.path.join(path2Xref, 'hg19_xref.txt'), sep='\t')
idx2gene = {int(idx): symbol for idx, symbol in zip( Xtable['EntrezGene ID'], Xtable['Approved Symbol']) if not math.isnan(idx)}

def go2genes( go, sep = ';' ):
    indices = go2idx.get(go)
    if indices is None:
        return ''
    else:
        return sep.join([idx2gene.get(idx) for idx in indices if idx2gene.get(idx) is not None])

prefix = "/Users/erickernfeld/Desktop/scRNA_data_analysis/cropseq_analysis/mimosca/results/after improved data cleaning/DE_TERA"
go_big_or_go_home = pd.read_csv(os.path.join(prefix, "GO.csv"), index_col = 0, header = 0)

try:
    os.makedirs(os.path.join(prefix, "GO_by_guide"))
except:
    pass

for guide in [re.sub("_Terarep.", "", guide_rep) for guide_rep in go_big_or_go_home.columns]:
    cols_gimme = [guide + "_Terarep1", guide + "_Terarep2"]
    cols_gimme = [colname for colname in cols_gimme if colname in go_big_or_go_home.columns]
    go_small = go_big_or_go_home.loc[:, cols_gimme].copy()
    go_small.loc[:, "name"]      = [ rowname.split('__')[1] for rowname  in go_small.index ]
    go_small.loc[:, "accession"] = [ rowname.split('__')[0] for rowname  in go_small.index ]
    go_small.loc[:, "genes"]     = [ go2genes( accession ) for accession in go_small.loc[:, "accession"] ]
    go_small.index = go_small.loc[:, "accession"]
    go_small = go_small.loc[:, ["accession", "name"] + cols_gimme + ["genes"]]
    pval_colnames = ["signed_log10_pval" + "_" + guide_rep[-4:] for guide_rep in cols_gimme]
    go_small.columns = ["accession", "name"] + pval_colnames + ["genes"]
    go_small = go_small.sort_values(pval_colnames[0], ascending = False)
    keep = [any( [abs(float(p)) >= 0.05 for p in pvals[1]] ) for pvals in go_small.loc[:, pval_colnames].iterrows() ]
    go_small = go_small.iloc[keep, :]
    go_small.to_csv(os.path.join( prefix, "GO_by_guide", guide + ".csv" ), index = False, quoting = True )