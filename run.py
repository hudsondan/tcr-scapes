import pandas as pd
import os, argparse
from benchmark import run
from functions.util_functions import get_time

def parse(args):
    '''Parse input arguments
    :param args: parameter inputs from terminal
    :return argdict: parsed arguments for benchmarking
    :rtype: dict'''

    # Confirm valid model selection

    if args.model_selection not in ['fast',
                                    'all',

                                    'clustcr',
                                    'hamming',
                                    'GIANA',
                                    'gliph2',
                                    'ismart',
                                    'length',
                                    'tcrdist',
                                    'tcrdist3',
                                    'vcluster',
                                    ]:
        raise AssertionError('Choose a model selection from "all","fast","clustcr","hamming","GIANA","gliph2","ismart","tcrdist","tcrdist3", "vcluster')
        
    elif args.model_selection=='all':
        model_selection = ['length',
                           'vcluster',
                            'clustcr',
                            'hamming',
                            'GIANA',
                            'gliph2',
                            'ismart',
                            'tcrdist3',
                            ]
    elif args.model_selection=='tcrdists':
        model_selection = ['tcrdist',
                            'tcrdist3',
                            ]
    else:
        model_selection = [args.model_selection]
        
    if args.chain_selection not in ['alpha','beta','both']:
        raise AssertionError('Choose a chain selection from ["alpha","beta", "both"]')
    
    if args.chain_selection in ['alpha','beta']:
       chain_selection = [args.chain_selection ]
    else:
        chain_selection=['alpha','beta']
        
    if args.paired in [None, 'None', False, 'False']:
        paired = None
    else:
        paired = True
    
    argdict  = {'Datetime': get_time(),
                'Input':args.input,
                'Experiment': args.expt,
                'CPUs': args.cpus,
                'Chain_selection':chain_selection,
                'Model_selection':model_selection,
                'Deduplicate':args.dedup,
                'Downsample': args.dsample,
                'Repeats':args.repeats,
                'N_olga':args.n_olga,
                'Min_clustsize': args.min_eps,
                'Paired':paired,
                'tcrdist_Method':args.tcrdist_method,
                'tcrdist_Radius':args.tcrdist_radius,
                'tcrdist_Hyper': args.tcrdist_hyper,
                'tcrdist_Subset': args.tcrdist_subset,
                'tcrdist_Chunk': args.tcrdist_chunk,
                'Save':args.save,
                }
    return argdict
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='ClustOx: Benchmarking Unsupervised Clustering Models for inference of TCR epitope specificity')
    parser.add_argument('-i', '--input', required=False, type=str, default='data/vdjdb.csv',
                        help='Set the input dataset')
    parser.add_argument('-c', '--cpus', required=False, type=int, default=1,
                        help='Set the number of available CPUs')
    parser.add_argument('-cs', '--chain_selection', required=False, type=str, default = 'beta',
                        help='Select chain selection from ["alpha","beta","both"]')
    parser.add_argument('-m', '--model_selection', required=False, type=str, default = 'all',
                        help='Select models from ["fast" (ClusTCR, GIANA, Hamming, Length, Random),"all" (ClusTCR, GIANA, GLIPH2, Hamming, iSMART, tcrdist3, Length, Random)]')
    parser.add_argument('-p', '--paired', required=False, type=bool, default = True,
                        help='Select whether to retain only paired chain TCRs')
    parser.add_argument('-ds', '--dsample', required=False, type=int, default = 1000,
                        help='If minimum clusters enabled, downsample to N TCRs per epitope')
    parser.add_argument('-d', '--dedup', required=False, type=bool, default = True,
                        help='Remove duplicates on a single chain selection')
    parser.add_argument('-r', '--repeats', required=False, type=int, default = 1,
                        help='Number of experimental repeats')
    parser.add_argument('-no', '--n_olga', required=False, type=int, default=0,
                        help='Number of randomly generated TCRs to spike in')
    parser.add_argument('-me', '--min_eps', required=False, type=int, default = 1000,
                        help='Retain epitopes having a minimum number of representatives')
    parser.add_argument('-tcdm', '--tcrdist_method', required=False, type=str, default = 'DBSCAN',
                        help='Choose tcrdist3 clustering method from ["DBSCAN","KMeans"]')
    parser.add_argument('-tcdr', '--tcrdist_radius', required=False, type=int, default = 50,
                        help='Choose tcrdist3 metaclonotype radius')
    parser.add_argument('-tcdh', '--tcrdist_hyper', required=False, type=float, default = 0.5,
                        help='Choose tcrdist3 clustering model hyperparameter (DBSCAN: n_eps, KMeans n_clusts)')
    parser.add_argument('-tcds', '--tcrdist_subset', required=False, type=bool, default = True,
                        help='Select only TCRs with a corresponding gene in the tcrdist lookup')
    parser.add_argument('-tcdc', '--tcrdist_chunk', required=False, type=bool, default = True,
                        help='Chunk tcrdist3 measurements to avoid OOM')
    parser.add_argument('-s', '--save', required=False, type=bool, default = None,
                        help='Save output clusters')
    parser.add_argument('-ex', '--expt', required=False, type=str, default = '',
                        help='Set experiment name')
    
    argdict = parse(parser.parse_args()) # Parse input arguments

    data = pd.read_csv(argdict['Input'],low_memory=False)    # Load data

    # Prepare output directories

    if argdict['Save']:
        if argdict['Datetime'] not in os.listdir('results'):
            os.mkdir("results/%s"%(argdict['Datetime']))

    # Run
    run(data,argdict)