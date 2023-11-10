
import time, os, sys
from tcrdist.repertoire import TCRrep
from sklearn.cluster import DBSCAN, KMeans
from functions.processing import get_bio
import pandas as pd
import numpy as np
import networkx as nx

import multiprocessing as mp
from functions.processing import get_bio
from clustcr import Clustering
from functions.tcrdist_cpp import *
# from functions.tcrdist_cpp import cluster_TCRDist_matrix_cpp

from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory2
# from tcrdist.rep_funcs import  compute_n_tally_out_of_memory2

def run_clustcr(df, chain_selection,cpus):
    '''Run ClusTCR clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :param cpus: # CPUs
    :type cpus: int
    :return mapper: cluster results
    :rtype mapper: dictionary
    :return t: runtime (s)
    :rtype t: float
    '''

    t0 = time.time()
    clustering = Clustering(n_cpus=cpus)
    mapper = {cdr3: cluster for cdr3, cluster in clustering.fit(df['cdr3.%s'%(chain_selection)], v_gene_col=df['v.%s'%(chain_selection)]).clusters_df.values.tolist()}
    t1 = time.time()
    return mapper, t1-t0

def run_GIANA(input_df, chain_selection, cpus):
    '''Run GIANA clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :param cpus: # CPUs
    :type cpus: int
    :return mapper: cluster results
    :rtype mapper: dictionary
    :return t: runtime (s)
    :rtype t: float
'''

    # Reformat input for GIANA
    cols=['%s.%s'%(x, chain_selection) for x in ['cdr3','v',
                                                ]]
    df = input_df[cols]
    df.columns=['CDR3','V']

    # Run GIANA
    wdir = os.getcwd()
    GIANA_PATH = os.path.join(wdir, 'functions/GIANA/')
    os.chdir(GIANA_PATH)
    df.to_csv('input.txt', index=False, header=False, sep='\t')
    print('Clustering {} sequences with GIANA.'.format(len(df)))

    # Perform GIANA algorithm on test sequences
    t0 = time.time()
    os.system('python GIANA4.1.py -f input.txt -O input_clustered.txt -v True -N {}'.format(cpus))
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    with open(os.path.join(GIANA_PATH, 'input_clustered.txt'), 'r') as f:
        clusters = f.read().splitlines()[3:]
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3',
                                                                            'cluster',
                                                                            'V', 
                                                                            ])

    # Map from outout to input sequences using bioidentities
    clusters = get_bio(clusters.rename(columns={'CDR3':'cdr3.%s'%(chain_selection),
                                                'V': 'v.%s'%(chain_selection),
                                                }),chain_selection,
                                                False,
                                                )
    clustmap = {bio: cluster for bio, cluster in clusters[['bio','cluster']].values.tolist()}

    os.chdir(wdir)

    return clustmap, t

def run_ismart(input_df, chain_selection,cpus):
    '''Run ISMART clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :param cpus: # CPUs
    :type cpus: int
    :return mapper: cluster results
    :rtype mapper: dictionary
    :return t: runtime (s)
    :rtype t: float
    '''

    # Reformat input for iSMART
    cols=['%s.%s'%(x, chain_selection) for x in ['cdr3','v',
                                                 ]]
    df = input_df[cols]
    df.columns=['CDR3',
                'V',
                ]

    # Run iSMART
    wdir = os.getcwd()
    ISMART_PATH = os.path.join(wdir, 'functions/ismart/lib')
    os.chdir(ISMART_PATH)
    df.to_csv('input.txt', index=False, header=False, sep='\t')
    print('Clustering {} sequences with iSMART.'.format(len(df)))
    t0 = time.time()
    os.system('python iSMARTf3.py -f input.txt -v True -N {}'.format(cpus))
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    with open(os.path.join(ISMART_PATH, 'input_clustered_v3.txt'), 'r') as f:
        clusters = f.read().splitlines()[3:]
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3', 
                                                                            'V', 
                                                                            'cluster'])

    # Map from outout to input sequences using bioidentities
    clusters = get_bio(clusters.rename(columns={'CDR3':'cdr3.%s'%(chain_selection),
                                                'V': 'v.%s'%(chain_selection),
                                                }),
                                                chain_selection,
                                                            False)

    clustmap = {bio: cluster for bio, cluster in clusters[['bio','cluster']].values.tolist()}
    os.chdir(wdir)

    return clustmap, t

def run_gliph2(df, chain_selection):
    '''Run GLIPH2 clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :return mapper: cluster results
    :rtype mapper: dictionary
    :return t: runtime (s)
    :rtype t: float
    ''' 

    # Prepare input for GLIPH

    cols = ['CDR3', 'V', 'J', 'CDR3_alpha',
            'subject:condition', 'count','epitope','bio']

    if 'subject:condition' not in df.columns:
        df['subject:condition'] = ['none:none']*len(df)
    if 'count' not in df.columns:
        df['count']=np.ones(len(df)).astype(int)

    if chain_selection == 'alpha':
        df = df.rename(columns={'cdr3.alpha': 'CDR3',
                                                    'v.alpha': 'V',
                                                    'j.alpha': 'J'})

    elif chain_selection == 'beta':
        df = df.rename(columns={'cdr3.beta': 'CDR3',
                                                'v.beta': 'V',
                                                'j.beta':'J',
                                                })

    else:
        df = df.rename(columns={'cdr3.beta': 'CDR3',
                                    'v.beta': 'V',
                                    'j.beta':'J',
                                    'cdr3.alpha':'CDR3_alpha',
                                                })
    colsx = [c for c in cols if c in df.columns]
    df=df[colsx]
    wdir= os.getcwd()

    # Run GLIPH2

    # NB: .centos executable required for Linux Centos, .OSX for Mac

    print('\n*** Clustering %s %s chains with GLIPH2 **'%(len(df),chain_selection))
    GLIPH2_PATH = os.path.join(wdir,'functions/gliph2/lib')
    os.chdir(GLIPH2_PATH)
    df.to_csv(os.path.join(GLIPH2_PATH,'metarepertoire.txt'), index=False, header=False, sep='\t')
    t0 = time.time()
    if sys.platform.lower() == 'darwin':
        os.system('./irtools.osx -c parameters.txt')
    elif sys.platform.lower() == 'linux':
        os.system('./irtools.centos -c parameters.txt')
    else:
        raise SystemError('GLIPH2 can run on Mac (Darwin) or Linux systems.',
                            'Windows users are directed to the GLIPH2 webtool, a link to which is provided in the readme',
                            'or to select a different model')
    t1 = time.time()
    t = t1 - t0

    # Reformat gliph2 clustering results

    clusters = pd.DataFrame()
    with open('metarepertoire_cluster.txt', 'r') as f:
        results = f.read().splitlines()
    c = 0
    for line in results:
        columns = line.split(' ')
        p = columns[0]
        motif = columns[3]
        cluster = columns[4:]
        if len(cluster) >= 2:
            nodes = pd.DataFrame({'CDR3': cluster})
            nodes['cluster'] = c
            nodes['motif'] = motif
            nodes['p'] = p
            clusters = pd.concat([clusters,nodes])
            c += 1

    os.chdir(wdir)
    mapper = {cdr3: clusters[clusters['CDR3']==cdr3]['cluster'].iloc[0] for cdr3 in clusters['CDR3'].unique()}

    return mapper, t

def hamming_hash(cdr3):
    '''Hamming distance computation with efficient hashing,
    :param cdr3: list of sequences
    :type cdr3: Series
    :return edgedict: dictionary of matched sequences
    :rtype edgedict: dict'''

    print('Calculating Hamming distances')

    cdr3hash = dict()
    for cdr in cdr3:
        for hash in (cdr[::2], cdr[1::2]):
            if hash not in cdr3hash:
                cdr3hash[hash] = set()
            cdr3hash[hash].add(cdr)
    
    edgedict = {}
    for i, hash in enumerate(cdr3hash):
        if len(cdr3hash[hash]) >= 1:
            for cdr1 in cdr3hash[hash]:
                if cdr1 not in edgedict.keys():
                    edgedict[cdr1] = i
                    for cdr2 in cdr3hash[hash]:
                        if cdr2 not in edgedict.keys():
                            if cdr1 != cdr2:
                                if cdr1 <= cdr2:
                                    if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) == 1:
                                        edgedict[cdr2] = i
    
    return edgedict


def run_hamming(df,chain_selection):
    '''Run Hamming Distance clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :return mapper: cluster results
    :rtype mapper: dictionary
    :return t: runtime (s)
    :rtype t: float
    ''' 

    cdr3=df['cdr3.%s'%(chain_selection)].unique()
    print('\n*** Clustering %s %s chains on Hamming distance **' % (len(df),
                                                                    chain_selection))

    t0 = time.time() # Initialise runtime
    mapper = hamming_hash(cdr3)
    t1 = time.time()
    t = t1 - t0
    
    return mapper, t

def shuffle(array):
    '''Randomly shuffle an input array'''

    idx = np.arange(len(array))
    np.random.shuffle(idx)
    out = array[idx]
    assert len(out)==len(array)
    assert not np.array_equal(out,array)
    return array[idx]

def run_random(df):
    '''Run random clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :return df: cluster results
    :rtype df: Pandas DataFrame
    :return t: runtime (s)
    :rtype t: float
    ''' 

    print('Random clustering')

    t0 =time.time()
    clustdict = {ep: id for id, ep in enumerate(df['epitope'].unique())}
    df.loc[:,'cluster']=shuffle(df['epitope'].values)
    df.loc[:,'cluster']=df['cluster'].map(clustdict)
    t1=time.time()

    return df, t1-t0

def run_vcluster(df, pairing):

    t0 = time.time()
    mapper = {v:i for i,v in enumerate(df['v.%s'%(pairing)].unique())}
    df['cluster']=df['v.%s'%(pairing)].map(mapper)
    t1 = time.time()

    return df, t1-t0





def cluster_TCRDist_matrix(S, seqs, method,cpus=1,hyperparam=None):
    '''Cluster distance matrix from tcrdist
    :param S: Distance matrix from tcrdist3
    :type S: array
    :param seqs: input data
    :type seqs: Pandas DataFrame
    :param method: clustering method
    :type method: str
    :param cpus: # CPUs
    :type cpus: int
    :param hyperparam: hyperparameter for clustering method
    :type hyperparam: float
    :return mapper: cluster results
    :rtype mapper: dictionary

    '''
    # Available methods
    methods = ['DBSCAN', 'KMEANS']
    method=method.upper()
    assert method in methods, r'Please choose one of the following: /n %s' % methods

    if method == 'DBSCAN':
        if not hyperparam:
            hyperparam=0.5
        clustering = DBSCAN(eps=hyperparam, min_samples=2, n_jobs=cpus).fit(S)
        labels = clustering.labels_
        
    elif method == 'KMEANS':
        os.environ["OMP_NUM_THREADS"] = "%s"%(str(cpus)) # export OMP_NUM_THREADS
        if not hyperparam:
            hyperparam=500
        kmeans = KMeans(init='random',
                        n_init=10,
                        n_clusters=int(hyperparam)).fit(S)
        labels = kmeans.labels_

    return {seq: label for seq, label in zip(seqs['bio'].values,labels) if label!=-1}

def run_tcrdist3(df, chain_selection, cpus, method = 'DBSCAN', radius=50, hyper=None, chunk=True):
    '''Run tcrdist3 clustering algorithm
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :param cpus: # CPUs
    :type cpus: int
    :param method: clustering method
    :type method: str
    :param radius: tcrdist metaclonotype radius
    :type radius: int
    :param hyper: clustering algorithm hyperparameter
    :type hyper: float
    :return df2: cluster results
    :rtype df2: pandas DataFrame
    :return t: runtime (s)
    :rtype t: float
    '''
    
    
    # Reformat input dataframe for tcrdist

    df2 = df.copy()
    cdr3a = 'cdr3_a_aa'
    va = 'v_a_gene'
    ja = 'j_a_gene'
    cdr3b = 'cdr3_b_aa'
    vb = 'v_b_gene'
    jb = 'j_b_gene'

    if not 'count' in df2.columns:
        df2['count']=[1]*len(df2)

    if chain_selection == 'alpha':
        df2 = df2.rename(columns={'cdr3.alpha': cdr3a,
                                                'v.alpha': va,
                                                'j.alpha':ja})
        df_epi = df2[[cdr3a, va, ja, 'bio','epitope', 'count']]

    elif chain_selection == 'beta':
        df2 = df2.rename(columns={'cdr3.beta': cdr3b,
                                                'v.beta': vb,
                                                'j.beta':jb})
        df_epi = df2[[cdr3b, vb, jb, 'bio','epitope', 'count']]

    else:
        df2 = df2.rename(columns={'cdr3.alpha': cdr3a,
                                                'v.alpha': va,
                                                'j.alpha':ja,
                                                'cdr3.beta': cdr3b,
                                                'v.beta': vb,
                                                'j.beta':jb})
        df_epi = df2[[va, cdr3a, ja, vb, cdr3b, jb, 'bio','epitope', 'count']]

    seqs = df_epi.drop(columns=['epitope'], axis=1).reset_index(drop=True)

    if chain_selection in ['alpha', 'beta']:
        chain = [chain_selection]
    else:
        chain = ['alpha', 'beta']

    # Run tcrdist3

    print('\n*** Clustering %s %s chains with tcrdist3' % (len(seqs), chain_selection))

    t0 = time.time()

    # Create tcrdist object
    tr = TCRrep(cell_df=seqs,   # Input data
                organism='human',   # Organism
                chains=chain,       # Chain selection
                infer_all_genes=True,   # Infer genes
                infer_cdrs=True,        # Infer CDRs
                infer_index_cols=True,  # Infer column names
                store_all_cdr=True,     # Store CDRs
                deduplicate=False,      # Deduplicate
                compute_distances=False)    # Compute distances
    
    # Compute tcrdist distances using sparse rect distances


    tr.cpus = cpus

    if chunk:
        if chain ==['alpha']:
            name = 'alpha'
        else:
            name = 'beta'

        S, _ = compute_pw_sparse_out_of_memory2(tr = tr,
            row_size      = 50,
            pm_processes  = cpus,
            pm_pbar       = True,
            max_distance  = radius,
            reassemble    = True,
            cleanup       = True,
            assign        = True)
        S=S[name]
    else:
        tr.compute_sparse_rect_distances(radius=radius, chunk_size=500)
        if chain_selection == 'alpha':
            S = tr.rw_alpha
        elif chain_selection == 'beta':
            S = tr.rw_beta
        else:
            S = tr.rw_beta

    # Convert columns back to input
    namedict = {cdr3a:'cdr3.alpha',
            va:'v.alpha',
            ja:'j.alpha',
    cdr3b: 'cdr3.beta',
    vb: 'v.beta',
    jb: 'j.beta'}
        
    # Record results
    t1 = time.time()
    print('Distance matrix calculated in %s seconds, clustering with %s' % (t1 - t0, method))
    t2 = time.time()

    # Cluster tcrdist matrix
    mapper = cluster_TCRDist_matrix(S, seqs, method = method, cpus=cpus, hyperparam=hyper)
    t3 = time.time()

    tdist = t1-t0
    tclust = t3 - t2

    # map cluster back to sequence based on bioidentity
    df2.loc[:,'cluster']=df2['bio'].map(mapper)
    df2=df2.rename(columns=namedict)

    return df2, [tdist,tclust]

def run_tcrdistcpp(input_df, chain_selection, cpus=1, method = 'DBSCAN',radius = 50, hyper = None):
    '''Run tcrdist C++ clustering algorithm
    Adapted from https://github.com/phbradley/conga
    :param input_df: input data
    :type input_df: Pandas DataFrame
    :param chain_selection: TCR chain selection (alpha/beta)
    :type chain_selection: str
    :param cpus: # CPUs
    :type cpus: int
    :param method: clustering method
    :type method: str
    :param radius: tcrdist metaclonotype radius
    :type radius: int
    :param hyper: clustering algorithm hyperparameter
    :type hyper: float
    :return input_df: cluster results
    :rtype input_df: pandas DataFrame
    :return t: runtime (s)
    :rtype t: float
    '''
    # Prepare data for tcrdist processing 
    if chain_selection =='alpha':
        tcrs = [((input_df.iloc[i]['v.alpha'],None, input_df.iloc[i]['cdr3.alpha']), (None,None,None)) for i in range(len(input_df))]
    elif chain_selection =='beta':
        tcrs = [((None,None, None), (input_df.iloc[i]['v.beta'],None, input_df.iloc[i]['cdr3.beta'])) for i in range(len(input_df))]
    else:
        tcrs = [((input_df.iloc[i]['v.alpha'],None, input_df.iloc[i]['cdr3.alpha']), (input_df.iloc[i]['v.beta'],None,input_df.iloc[i]['cdr3.beta'])) for i in range(len(input_df))]

    input_df.loc[:,'bio']=tcrs
    uniques = list(set(tcrs))
    t0 = time.time()
    
    # Compute tcrdist distance matrix
    D2 = calc_tcrdist_matrix_cpp(uniques,'human',
                                 single_chain=chain_selection,
                                 threshold = radius)

    t1=time.time()
    print('Distance matrix calculated in %s seconds, clustering with %s' % (t1 - t0, method))
    t2=time.time()

    # Cluster tcrdist distance matrix
    mapper = cluster_TCRDist_matrix_cpp(D2, 
                                    uniques,
                                    method,
                                    cpus=cpus,
                                    hyperparam=hyper)
    t3 = time.time()

    # Map from sequences to cluster allocation
    input_df.loc[:,'cluster']=input_df['bio'].map(mapper)

    tdist = t1-t0
    tclust = t3 - t2
    return input_df, [tdist,tclust]
