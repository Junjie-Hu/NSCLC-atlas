#!/usr/bin/env python
import time
import resource
import os, sys, pandas, numpy, pathlib, getopt
import CytoSig

from scipy import stats
from statsmodels.stats.multitest import multipletests
from optparse import OptionParser

usage = "Usage:\nTres.py -i <input file, default: Treg_a_sample_CytoSig.csv> -n <normalize, default: False> -o <output,default: Treg_a_proliferation.csv>"
optparser = OptionParser(usage=usage)
optparser.add_option("-i", "--input", help="file to process")
optparser.add_option("-n", "--normalize", help="normalize or not (True or False)")
optparser.add_option("-o", "--output", help="output name")
(options, args) = optparser.parse_args(sys.argv)

if not options.input:
	optparser.print_help()
	sys.exit(-1)

def profile_geneset_signature(expression):
    signature = []
    fin = open(os.path.join(fpath, 'Tres.kegg'))
    for l in fin:
        fields = l.strip().split('\t')
        
        s = fields[2:]
        signature.append(pandas.Series(numpy.ones(len(s)), index=s, name=fields[0]))
    fin.close()
    
    signature = pandas.concat(signature, axis=1, join='outer', sort=False)
    signature.fillna(0, inplace=True)
    
    common = expression.index.intersection(signature.index)
    signature, expression = signature.loc[common], expression.loc[common]
    
    background = signature.mean(axis=1)
    background.name = 'study bias'
    
    X = signature.loc[:, ['KEGG_CELL_CYCLE', 'KEGG_DNA_REPLICATION']].mean(axis=1)
    X.name = 'Proliferation'
    
    X = pandas.concat([background, X], axis=1, join='inner')
    
    result = CytoSig.ridge_significance_test(X, expression, alpha=0, alternative="two-sided", nrand=0, verbose=verbose_flag)
    
    return result[2].loc['Proliferation']

err_tol = 1e-8
verbose_flag = False
fpath = "/media/inspur/AS2150G2/HJJ/scrna/CytoSig/Tres/src"

input_file = options.input
normalization = options.normalize

expression = pandas.read_csv(input_file, sep=',', index_col=0)

# normalization
if normalization:
    background_expression = expression.mean(axis=1)
    expression = expression.subtract(background_expression, axis=0)

result_prolifertion = profile_geneset_signature(expression)

## save 
result_prolifertion.to_csv(options.output)



