import pandas as pd
import multiprocessing as mp
import os, time, argparse

from functions.util_functions import get_time
from functions.processing import *
from functions.models import run_clustcr, run_hamming,  run_GIANA, run_gliph2, run_ismart, run_random, run_tcrdist3, run_tcrdistcpp
from functions.metrics import *

class Benchmark:
    
    def __init__(self, argdict):
        '''Initialise benchmark cluster class
        :param argdict: dictionary of input parameters
        :type argdict: dictionary'''
        
        self.argdict = argdict
        self.cpus = argdict['CPUs']
        self.chain = argdict['Chain']
        self.datetime = argdict['Datetime']
        self.model = argdict['Model']
        self.save = argdict['Save']
        
    def baseline(self, data, 
                clusterscores, 
                epscores, 
                statistics):
        '''Establish a random baseline for comparison of model performance
                :param data: source data
        :type data: DataFrame
        :param clusterscores: set of global cluster scores
        :type clusterclores: DataFrame
        :param epscores: set of epitope-specific cluster scores
        :type epscores: DataFrame
        :param statistics: cluster statistics
        :type statistics: DataFrame
        '''
    
        self.argdict['N_total'] = len(data)
        data, randtime = run_random(data)   
        self.argdict['Model'] = 'random'     
        self.argdict['Runtime'] =randtime
        self.argdict['N_clusters']= len(data['cluster'].dropna().unique())
        c_scores, e_scores, stats = score(data,self.argdict)
        
        clusterscores=pd.concat([clusterscores,c_scores])
        epscores=pd.concat([epscores,e_scores])
        statistics=pd.concat([statistics,stats])
        
        if self.save:
            data.to_csv('results/%s/%s_%s.csv'%(self.datetime,'random',self.chain))

        return clusterscores, epscores, statistics
    
    def getscore(self, data,
                    clusterscores,
                    epscores,
                    statistics):
        
        '''Return evaluation of cluster and predictive metrics for a given model
        and parameter set
        :param data: source data
        :type data: DataFrame
        :param clusterscores: set of global cluster scores
        :type clusterclores: DataFrame
        :param epscores: set of epitope-specific cluster scores
        :type epscores: DataFrame
        :param statistics: cluster statistics
        :type statistics: DataFrame
        '''
        data2=data.copy()
        if self.model =='clustcr':
            mapper, t = run_clustcr(data2, 
                                    self.chain, 
                                    self.cpus)
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map(mapper)

        if self.model =='GIANA':
            data2 = get_bio(data2, 
                           self.chain, 
                           False)
            mapper, t = run_GIANA(data2,
                                  self.chain,
                                  cpus=self.cpus)
            data2.loc[:,'cluster']=data2['bio'].map(mapper)

        elif self.model =='gliph2':
            mapper, t = run_gliph2(data2, self.chain)
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map(mapper)
        
        elif self.model =='ismart':
            data2 = get_bio(data2, 
                           self.chain, 
                           False)
            mapper, t = run_ismart(data2, 
                                   self.chain,
                                   cpus=self.cpus)
            data2.loc[:,'cluster']=data2['bio'].map(mapper)

        elif self.model=='hamming':
            clusters, t = run_hamming(data2, self.chain)
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map(clusters)
            # cnts = [l for l,n in data2['cluster'].value_counts().reset_index()[['cluster','count']].values.tolist() if n<=1]
            # data2.loc[:,'cluster']=data2.loc[:,'cluster'].replace(cnts,np.nan)
        
        elif self.model=='length':
            print('Clustering on cdr3 length')
            t0 =time.time()
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map({cdr3: 0 if type(cdr3)==float else len(cdr3) for cdr3 in data2['cdr3.%s'%(self.chain)].values})
            # cnts = [l for l,n in data2['cluster'].value_counts().reset_index()[['cluster','count']].values.tolist() if n<=1]
            # data2.loc[:,'cluster']=data2.loc[:,'cluster'].replace(cnts,np.nan)
            t1 = time.time()
            t=t1-t0
            
        elif self.model =='tcrdist':

            data2, t = run_tcrdistcpp(data2,
                                    self.chain,
                                    cpus=self.cpus,
                                    method=self.argdict['tcrdist_Method'],
                                    radius=self.argdict['tcrdist_Radius'],
                                    hyper=self.argdict['tcrdist_Hyper'],
                                    )
            
            # data2.loc[:,'cluster']=data2['bio'].map({bio:cluster for (bio, cluster) in clusters[['bio','cluster']].values.tolist()})
            
        elif self.model =='tcrdist3':
            data2 = get_bio(data2, 
                           self.chain, 
                           True)
            data2, t = run_tcrdist3(data2,
                                    self.chain,
                                    cpus=self.cpus,
                                    method=self.argdict['tcrdist_Method'],
                                    radius=self.argdict['tcrdist_Radius'],
                                    hyper=self.argdict['tcrdist_Hyper'],
                                    )
            
            # data2.loc[:,'cluster']=data2['bio'].map({bio:cluster for (bio, cluster) in clusters[['bio','cluster']].values.tolist()})
        # Drop singlets
        cnts = [l for l,n in data2['cluster'].value_counts().reset_index()[['cluster','count']].values.tolist() if n<=1]
        data2.loc[:,'cluster']=data2.loc[:,'cluster'].replace(cnts,np.nan)
        
        if self.save:
            data2.to_csv('results/%s/%s_%s.csv'%(self.datetime,self.model,self.chain))
            
        self.argdict['Runtime'] =t
        self.argdict['N_clusters']= len(data2['cluster'].dropna().unique())
        c_scores, e_scores, stats = score(data2, self.argdict)
        
        clusterscores=pd.concat([clusterscores,c_scores])
        epscores=pd.concat([epscores,e_scores])
        statistics=pd.concat([statistics,stats])
        
        clusterscores, epscores, statistics = self.baseline(data.copy(),
                                                            clusterscores,
                                                            epscores,
                                                            statistics)
        return clusterscores, epscores, statistics
    
def run(input_data, argdict):
    c, e, s = [pd.DataFrame(), pd.DataFrame(),pd.DataFrame()]
    minclusts=argdict['Min_clustsize']
    datetime = argdict['Datetime']
    
    if argdict['Deduplicate']:
        data=input_data.drop_duplicates(subset=['cdr3.alpha',
                                         'v.alpha',
                                         'j.alpha',
                                         'cdr3.beta',
                                         'v.beta',
                                         'j.beta'])
        
    for dataset in data['dataset'].unique():
        
        sub=data[data['dataset']==dataset].copy()

        sub = generate_data(sub,
                            chain_selection=argdict['Chain_selection'],
                            n=minclusts,
                            paired = argdict['Paired'],
                            tcrd_sub=argdict['tcrdist_Subset'])
        sub.to_csv('results/input_data_vdj.csv')
        argdict['Dataset']=dataset
        for i in range(argdict['Repeats']):
            if argdict['Downsample']>0:
                
                df = pd.concat([sub[sub['epitope']==ep].sample(frac=1).iloc[:argdict['Downsample']] for ep in sub['epitope'].unique()])
            else:
                df=sub.copy()
            for chain_selection in argdict['Chain_selection']:
                df2=df.copy()

                argdict['Chain']= chain_selection
                if argdict['N_olga'] > 0:
                    olga = prepare(get_olga(chain_selection,
                                            argdict['N_olga']),
                                   chain_selection=argdict['Chain_selection'],
                                   tcrd_sub = True)
                else:
                    olga=None
            
                if argdict['N_olga'] > 0:
                    df2=pd.concat([df2,olga])
                
                for model in argdict['Model_selection']:
                    df3=df2.copy()
                    argdict['N_total']= len(df3)
                    argdict['Model']=model
                    benchmark = Benchmark(argdict)
                    c, e, s = benchmark.getscore(df3,
                                        c,
                                        e,
                                        s)
    
    if not datetime in os.listdir('results'):
        os.mkdir('results/%s'%(datetime))
    c.to_csv('results/%s/total.csv'%(datetime))
    e.to_csv('results/%s/eps.csv'%(datetime))
    s.to_csv('results/%s/stats.csv'%(datetime))

def parse(args):
    if args.model_selection not in ['fast',
                                    'all',
                                    'tcrdists',
                                    'length',
                                    'clustcr',
                                    'hamming',
                                    'GIANA',
                                    'gliph2',
                                    'ismart',
                                    'tcrdist',
                                    'tcrdist3',
                                    ]:
        raise AssertionError('Choose a model selection from "fast","all","tcrdists","length","clustcr","hamming","GIANA","gliph2","ismart","tcrdist3","tcrdist"')
    
    if args.model_selection =='fast':
        model_selection = ['length',
                            'clustcr',
                            'GIANA',
                            'Hamming'
                            'tcrdist'
                            ]
        
    elif args.model_selection=='all':
        model_selection = ['length',
                            'clustcr',
                            'hamming',
                            'GIANA',
                            'gliph2',
                            'ismart',
                            # 'tcrdist',
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
    
    argdict  = {'Datetime': get_time(),
                'Experiment': args.expt,
                'CPUs': args.cpus,
                'Chain_selection':chain_selection,
                'Model_selection':model_selection,
                'Deduplicate':args.dedup,
                'Downsample': args.dsample,
                'Repeats':args.repeats,
                'N_olga':args.n_olga,
                'Min_clustsize': args.min_eps,
                'Paired':args.paired,
                'tcrdist_Method':args.tcrdist_method,
                'tcrdist_Radius':args.tcrdist_radius,
                'tcrdist_Hyper': args.tcrdist_hyper,
                'tcrdist_Subset': args.tcrdist_subset,
                'Save':args.save,
                }
    return argdict
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='ClustOx: Benchmarking Unsupervised Clustering Models for inference of TCR epitope specificity')
    parser.add_argument('-c', '--cpus', required=False, type=int, default=1,
                        help='Set the number of available CPUs')
    parser.add_argument('-cs', '--chain_selection', required=False, type=str, default = 'beta',
                        help='Select chain selection from ["alpha","beta","both"]')
    parser.add_argument('-m', '--model_selection', required=False, type=str, default = 'all',
                        help='Select models from ["fast" (ClusTCR, GIANA, Hamming, Length, Random),"all" (ClusTCR, GIANA, GLIPH2, Hamming, iSMART, tcrdist3, Length, Random)]')
    parser.add_argument('-p', '--paired', required=False, type=bool, default = True,
                        help='Select whether to retain only paired chain TCRs')
    parser.add_argument('-ds', '--dsample', required=False, type=int, default = 500,
                        help='If minimum clusters enabled, downsample to N TCRs per epitope')
    parser.add_argument('-d', '--dedup', required=False, type=bool, default = True,
                        help='Remove duplicates on a single chain selection')
    parser.add_argument('-r', '--repeats', required=False, type=int, default = 1,
                        help='Number of experimental repeats')
    parser.add_argument('-no', '--n_olga', required=False, type=int, default=0,
                        help='Number of randomly generated TCRs to spike in')
    parser.add_argument('-me', '--min_eps', required=False, type=int, default = 500,
                        help='Retain epitopes having a minimum number of representatives')
    parser.add_argument('-tcdm', '--tcrdist_method', required=False, type=str, default = 'DBSCAN',
                        help='Choose tcrdist3 clustering method from ["DBSCAN","KMeans"]')
    parser.add_argument('-tcdr', '--tcrdist_radius', required=False, type=int, default = 50,
                        help='Choose tcrdist3 metaclonotype radius')
    parser.add_argument('-tcdh', '--tcrdist_hyper', required=False, type=float, default = 0.5,
                        help='Choose tcrdist3 clustering model hyperparameter (DBSCAN: n_eps, KMeans n_clusts)')
    parser.add_argument('-tcds', '--tcrdist_subset', required=False, type=bool, default = True,
                        help='Select only TCRs with a corresponding gene in the tcrdist lookup')
    parser.add_argument('-s', '--save', required=False, type=bool, default = None,
                        help='Save output clusters')
    parser.add_argument('-ex', '--expt', required=False, type=str, default = '',
                        help='Set experiment name')
    
    argdict = parse(parser.parse_args())

    root= 'data/'

    # data = pd.read_csv(os.path.join(root,'combined_deduplicated_clean_ags.csv'), low_memory=True, index_col=0).rename(columns={'antigen.epitope_clean':'epitope'})
        # mcpas = data[data['dataset']=='McPas-TCR']
    # mcpas.to_csv(os.path.join(root,'mcpas.csv'))
    # data = pd.read_csv(os.path.join(root,'mcpas.csv'), low_memory=False, index_col=0)
    # data = pd.read_csv(os.path.join(root,'combined_deduplicated_clean_ags.csv'), low_memory=True, index_col=0).rename(columns={'antigen.epitope_clean':'epitope'})
    # data=data[data['dataset'].isin(['VDJdb','10X'])]
    # data.to_csv(os.path.join(root, 'vdjdb_10x.csv'))
    # iedb = data[data['dataset']=='IEDB']
    # iedb.to_csv(os.path.join(root,'iedb.csv'))  

    # mira = data[data['dataset']=='MIRA']
    # mira.to_csv(os.path.join(root,'mira.csv'))

    # data = pd.read_csv(os.path.join(root,'vdjdb_10x.csv'), low_memory=True, index_col=0)
    # vdjdb = data[data['dataset']=='VDJdb']
    # vdjdb.to_csv(os.path.join(root,'vdjdb.csv'))  
    # print(data['epitope'].value_counts())
    
    data = pd.read_csv(os.path.join(root,'vdjdb.csv'), low_memory=False, index_col=0)
    # data = pd.read_csv(os.path.join(root,'mcpas.csv'), low_memory=True, index_col=0)
    # data = pd.read_csv(os.path.join(root,'mira.csv'), low_memory=True, index_col=0)
    # data = pd.read_csv(os.path.join(root,'iedb.csv'), low_memory=True, index_col=0)

    if argdict['Save']:
        if argdict['Datetime'] not in os.listdir('results'):
            os.mkdir("results/%s"%(argdict['Datetime']))
            
    run(data,argdict)