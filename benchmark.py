import pandas as pd
import os, time
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
                # epscores, 
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
        data, randtime = run_random(data)   # Run random baseline   
        self.argdict['Model'] = 'random'     
        self.argdict['Runtime'] =randtime
        self.argdict['N_clusters']= len(data['cluster'].dropna().unique())  # Retain only clustered instances
        c_scores, stats = score(data,self.argdict)    # Get scores
        # Add scores to output DataFrames
        clusterscores=pd.concat([clusterscores,c_scores])   
        statistics=pd.concat([statistics,stats])
        
        if self.save:
            # Record cluster output
            data.to_csv('results/%s/%s_%s.csv'%(self.datetime,'random',self.chain))

        return clusterscores, statistics
    
    def getscore(self, data,
                    clusterscores,
                    # epscores,
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

        # Run models, and map cluster outputs to DataFrame

        if self.model =='clustcr':
            # Run ClusTCR
            mapper, t = run_clustcr(data2, 
                                    self.chain, 
                                    self.cpus)
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map(mapper)

        if self.model =='GIANA':
            # Run GIANA
            data2 = get_bio(data2, 
                           self.chain, 
                           False)
            mapper, t = run_GIANA(data2,
                                  self.chain,
                                  cpus=self.cpus)
            data2.loc[:,'cluster']=data2['bio'].map(mapper)

        elif self.model =='gliph2':
            # Run GLIPH2
            mapper, t = run_gliph2(data2, self.chain)
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map(mapper)
        
        elif self.model =='ismart':
            # Run iSMART
            data2 = get_bio(data2, 
                           self.chain, 
                           False)
            mapper, t = run_ismart(data2, 
                                   self.chain,
                                   cpus=self.cpus)
            data2.loc[:,'cluster']=data2['bio'].map(mapper)

        elif self.model=='hamming':
            # Run Hamming
            clusters, t = run_hamming(data2, self.chain)
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map(clusters)
        
        elif self.model=='length':
            # Clust on cdr3 length
            print('Clustering on cdr3 length')
            t0 =time.time()
            data2.loc[:,'cluster']=data2['cdr3.%s'%(self.chain)].map({cdr3: 0 if type(cdr3)==float else len(cdr3) for cdr3 in data2['cdr3.%s'%(self.chain)].values})
            t1 = time.time()
            t=t1-t0
            
        elif self.model =='tcrdist':
            # Run C++ implementation of tcrdist from CoNGA package
            data2, t = run_tcrdistcpp(data2,
                                    self.chain,
                                    cpus=self.cpus,
                                    method=self.argdict['tcrdist_Method'],
                                    radius=self.argdict['tcrdist_Radius'],
                                    hyper=self.argdict['tcrdist_Hyper'],
                                    )
            
        elif self.model =='tcrdist3':
            # Run python implmentation of tcrdist3
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
            
        # Drop singlets
        cnts = [l for l,n in data2['cluster'].value_counts().reset_index()[['cluster','count']].values.tolist() if n<=1]
        data2.loc[:,'cluster']=data2.loc[:,'cluster'].replace(cnts,np.nan)
        
        # Record clusters
        if self.save:
            data2.to_csv('results/%s/%s_%s.csv'%(self.datetime,self.model,self.chain))

        # Record performance for model selected
        self.argdict['Runtime'] =t
        self.argdict['N_clusters']= len(data2['cluster'].dropna().unique())
        c_scores, stats = score(data2, self.argdict)
        clusterscores=pd.concat([clusterscores,c_scores])
        statistics=pd.concat([statistics, stats])
        
        # Generate random baseline results from same dataset
        clusterscores, statistics = self.baseline(data.copy(),
                                                            clusterscores,
                                                            # epscores,
                                                            statistics)
        return clusterscores, statistics


def run(input_data, argdict):
    '''Run benchmarking analysis
    :param input_data: input TCR instances
    :type input_data: DataFrame
    :param argdict: input parameters
    :type argdict: dict
    '''

    # Initialise results record
    c, e, s = [pd.DataFrame(), pd.DataFrame(),pd.DataFrame()]
    minclusts=argdict['Min_clustsize']
    datetime = argdict['Datetime']
    
    # Drop duplicates on CDR3A-TRAV-TRAJ-CDR3B-TRBV-TRBJ bioidentities
    if argdict['Deduplicate']:

        data=input_data.drop_duplicates(subset=['cdr3.alpha',
                                         'v.alpha',
                                         'j.alpha',
                                         'cdr3.beta',
                                         'v.beta',
                                         'j.beta'])
        
    for dataset in data['dataset'].unique():
        
        # Split datasets
        sub=data[data['dataset']==dataset].copy()

        # Preprocess input_data
        sub = preprocess(sub,
                            chain_selection=argdict['Chain_selection'],
                            n=minclusts,
                            paired = argdict['Paired'],
                            tcrd_sub=argdict['tcrdist_Subset'])
        argdict['Dataset']=dataset

        for i in range(argdict['Repeats']):

            if argdict['Downsample']>0:
                # Downsample to epitopes having a certain number of representatives
                df = pd.concat([sub[sub['epitope']==ep].sample(frac=1).iloc[:argdict['Downsample']] for ep in sub['epitope'].unique()])
            else:
                df=sub.copy()

            for chain_selection in argdict['Chain_selection']:
                df2=df.copy()
                argdict['Chain']= chain_selection

                if argdict['N_olga'] > 0:
                    # Add decoy sequences with OLGA (Sethna et al 2019)
                    olga = prepare(get_olga(chain_selection,
                                            argdict['N_olga']),
                                   chain_selection=argdict['Chain_selection'],
                                   tcrd_sub = True)
                    df2=pd.concat([df2,olga])

                else:
                    olga=None
                
                for model in argdict['Model_selection']:
                    
                    df3=df2.copy()
                    argdict['N_total']= len(df3)
                    argdict['Model']=model
                    benchmark = Benchmark(argdict) # Initialise
    
                    # Run benchmarking for models in the input selection
                    c, s = benchmark.getscore(df3,
                                        c,
                                        s)
    # Initialise results folder
    if not datetime in os.listdir('results'):
        os.mkdir('results/%s'%(datetime))

    # Save outputs
    c.to_csv('results/%s/total.csv'%(datetime)) # Performance values
    s.to_csv('results/%s/stats.csv'%(datetime)) # Cluster statistics

