
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

import time

def run_clustcr(df, chain_selection,cpus):
    t0 = time.time()
    clustering = Clustering(n_cpus=cpus)
    mapper = {cdr3: cluster for cdr3, cluster in clustering.fit(df['cdr3.%s'%(chain_selection)], v_gene_col=df['v.%s'%(chain_selection)]).clusters_df.values.tolist()}
    t1 = time.time()
    return mapper, t1-t0

def run_GIANA(input_df, chain_selection,cpus):
    """Run GIANA clustering algorithm
    :param cluster_object: Cluster object
    :type cluster_object: Object of class Cluster
    :return output: dictionary of cluster results
    :return t: runtime
    :return clusters: DataFrame of output clusters"""

    # Reformat input for GIANA
    cols=['%s.%s'%(x, chain_selection) for x in ['cdr3','v',
                                                #  'j'
                                                ]]
    df = input_df[cols]
    # df.columns=['CDR3','V','J']
    df.columns=['CDR3','V']

    # Run GIANA
    # print('\n*** Clustering %s %s chains with GIANA **' % (len(df),chain_selection))
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
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3', 'cluster','V' 
                                                                            # 'J',
                                                                            ])

    clusters = get_bio(clusters.rename(columns={'CDR3':'cdr3.%s'%(chain_selection),
                                                'V': 'v.%s'%(chain_selection),
                                                # 'J': 'j.%s'%(chain_selection),
                                                }),chain_selection,
                                                False,
                                                # True,
                                                )
    clustmap = {bio: cluster for bio, cluster in clusters[['bio','cluster']].values.tolist()}
    # output = ClusteringResult(join_cdr3_v(clusters)).metrics(join_cdr3_v(cluster_object.epitopes,
    #                 data_type='epitope_data')).summary()

    os.chdir(wdir)

    return clustmap, t

def run_ismart(input_df, chain_selection,cpus):
    """Run GIANA clustering algorithm
    :param cluster_object: Cluster object
    :type cluster_object: Object of class Cluster
    :return output: dictionary of cluster results
    :return t: runtime
    :return clusters: DataFrame of output clusters"""

    # Reformat input for iSMART
    cols=['%s.%s'%(x, chain_selection) for x in ['cdr3','v',
                                                #  'j',
                                                 ]]
    df = input_df[cols]
    df.columns=['CDR3',
                'V',
                # 'J',
                ]

    wdir = os.getcwd()
    ISMART_PATH = os.path.join(wdir, 'functions/ismart/lib')

    os.chdir(ISMART_PATH)
    df.to_csv('input.txt', index=False, header=False, sep='\t')

    print('Clustering {} sequences with iSMART.'.format(len(df)))

    # Perform iSMART algorithm on test sequences
    t0 = time.time()
    os.system('python iSMARTf3.py -f input.txt -v True -N {}'.format(cpus))
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    with open(os.path.join(ISMART_PATH, 'input_clustered_v3.txt'), 'r') as f:
        clusters = f.read().splitlines()[3:]
        clusters = pd.DataFrame([x.split('\t') for x in clusters], columns=['CDR3', 'V', 
                                                                            # 'J',
                                                                            'cluster'])

    clusters = get_bio(clusters.rename(columns={'CDR3':'cdr3.%s'%(chain_selection),
                                                'V': 'v.%s'%(chain_selection),
                                                # 'J': 'j.%s'%(chain_selection),
                                                }),chain_selection,
                                                        #    True,
                                                            False)

    clustmap = {bio: cluster for bio, cluster in clusters[['bio','cluster']].values.tolist()}
    os.chdir(wdir)

    return clustmap, t

def run_gliph2(df, chain_selection):

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
            
        # df_epi = ddf[['CDR3', 'V', 'J', 'subject:condition','count','epitope','bio']]

    elif chain_selection == 'beta':
        df = df.rename(columns={'cdr3.beta': 'CDR3',
                                                'v.beta': 'V',
                                                'j.beta':'J',
                                                })
        
        # df_epi = ddf[['CDR3', 'V', 'J', 'subject:condition','count','epitope','bio']]

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
    """Hamming distance computation with efficient hashing,
        adapted from clusTCR edge list production function
    :param cdr3: list of sequences
    :type cdr3: Series
    :return edgedict: dictionary of matched sequences
    :rtype edgedict: dict"""

    print('Calculating Hamming distances')
    # Adapted from clusTCR edge list production function

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
    """Run ClusTCR clustering algorithm
    :param cluster_object: Cluster object
    :type cluster_object: Object of class Cluster
    :return output: dictionary of cluster results
    :return t: runtime
    :return clusters: DataFrame of output clusters"""


    
    # if chain_selection !='paired':
    cdr3=df['cdr3.%s'%(chain_selection)].unique()
    print('\n*** Clustering %s %s chains on Hamming distance **' % (len(df),
                                                                    chain_selection))
    # Run comparison


    t0 = time.time() # Initialise runtime
    mapper = hamming_hash(cdr3)

    t1 = time.time()
    t = t1 - t0
    
    return mapper, t


def normalize(edge):
    '''
    Introduce normalization property on an edge that is represented as a tuple.
    The normalization property will constraint the ordering of two nodes that
    make up an edge. This prevents duplicated edges.

    Parameters
    ----------
    edge : tuple
        Tuple of length two, containing two nodes, represented as integers.

    Returns
    -------
    (n1, n2) : tuple
        Normalized edge.
        
    '''
    n1, n2 = edge
    if n1 > n2:
        n1, n2 = n2, n1
    return (n1, n2)


def greedy_clustering(dm, threshold):
    '''
    Greedy graph clustering approach that uses a fixed distance-threshold to 
    assign nodes to cluster. The algorithm starts by computing all pairs 
    of sequences that satisfy the predefined distance threshold (edge list). 
    Next, it finds the sequence with the highest degree (i.e. the most neighbors), 
    assigns this as the cluster centre, and removes it and its neighbors 
    from the edge list. This process is repeated until all sequences are clustered.

    Parameters
    ----------
    dm : numpy.array
        Distance matrix.
    threshold : int
        Distance threshold for defining network edges.

    Returns
    -------
    res : pandas.DataFrame
        Dataframe containing clustering results.

    '''

    edges = np.argwhere(dm<=threshold)
    # print(len(edges))
    edges = set(map(normalize, edges)) # Remove duplicated edges
    edges = np.array(list(edges)) # Edgelist to array
    # print(len(edges))
    
    cid = 0
    res = pd.DataFrame()
    
    while len(edges) > 0:
        
        G = nx.from_edgelist(edges)
        degrees = pd.DataFrame(G.degree(), columns=['node', 'degree'])
        degrees = degrees.set_index('node')
        degrees = degrees.sort_values(by='degree', ascending=False)
        max_degree = degrees.idxmax().values
        
        cluster = edges[np.where(edges[:,0]==max_degree)]
        ids = np.unique(cluster)
        cids = [cid] * len(ids)

        if len(ids) <= 1:
            break
        
        cluster_iter = pd.DataFrame({'seq_id':ids,'cluster':cids})
        res = pd.concat([res,cluster_iter])
        
        for i in ids:
            edges = edges[np.where(edges[:,0]!=i)] # Remove from column 1
            edges = edges[np.where(edges[:,1]!=i)] # Remove from column 2
        
        cid += 1
            
    return res

def shuffle(array):
    idx = np.arange(len(array))
    np.random.shuffle(idx)
    out = array[idx]
    assert len(out)==len(array)
    assert not np.array_equal(out,array)
    return array[idx]

def run_random(df):
    print('Random clustering')

    t0 =time.time()
    clustdict = {ep: id for id, ep in enumerate(df['epitope'].unique())}
    df.loc[:,'cluster']=shuffle(df['epitope'].values)
    df.loc[:,'cluster']=df['cluster'].map(clustdict)
    t1=time.time()
    return df, t1-t0

def cluster_TCRDist_matrix(S, seqs, method,cpus=1,hyperparam=None):
    '''
        Function for clustering of the TCRDist-calculated distance matrix.
        The function provides several methods for clustering, which include:
            - Greedy: fixed-distance threshold clustering method that groups
            of sequences that are connected in a graph.
            - DBSCAN: density-based clustering method that groups of densely
            packed points.
            - Kmeans: iterative clustering approach that partitions the data
            into k clusters.

        Parameters
        ----------
        dm : numpy.array, optional
            TCRDist-calculated distance matrix. The default is None.
        cdr3 : pandas.Series, optional
            pandas.Series containing the input CDR3 sequences. The default is None.
        gt : pandas.DataFrame, optional
            Ground truth data. The default is None.
        method : str, optional
            Clustering method used to cluster the TCRDist distance matrix. 
            Available methods include: greedy, DBSCAN and Kmeans.
            The default is 'DBSCAN'.

        Returns
        -------
        Clustering results

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
        
    # cnts = Counter(labels)
    # shortlist = [l for l in cnts.keys() if cnts[l]>=2]

    return {seq: label for seq, label in zip(seqs['bio'].values,labels) if label!=-1}

def run_tcrdist3(input_df, chain_selection, cpus, method = 'DBSCAN', radius=50, hyper=None):

    df = input_df.copy()

    cdr3a = 'cdr3_a_aa'
    va = 'v_a_gene'
    ja = 'j_a_gene'
    cdr3b = 'cdr3_b_aa'
    vb = 'v_b_gene'
    jb = 'j_b_gene'

    if not 'count' in df.columns:
        df['count']=[1]*len(df)

    if chain_selection == 'alpha':
        df = df.rename(columns={'cdr3.alpha': cdr3a,
                                                'v.alpha': va,
                                                'j.alpha':ja})
        
        df_epi = df[[cdr3a, va, ja, 'bio','epitope', 'count']]


    elif chain_selection == 'beta':
        df = df.rename(columns={'cdr3.beta': cdr3b,
                                                'v.beta': vb,
                                                'j.beta':jb})
        
        df_epi = df[[cdr3b, vb, jb, 'bio','epitope', 'count']]

    else:
        df = df.rename(columns={'cdr3.alpha': cdr3a,
                                                'v.alpha': va,
                                                'j.alpha':ja,
                                                'cdr3.beta': cdr3b,
                                                'v.beta': vb,
                                                'j.beta':jb})
        
        df_epi = df[[va, cdr3a, ja, vb, cdr3b, jb, 'bio','epitope', 'count']]

    seqs = df_epi.drop(columns=['epitope'], axis=1).reset_index(drop=True)

    if chain_selection in ['alpha', 'beta']:
        chain = [chain_selection]
    else:
        chain = ['alpha', 'beta']

    # Run tcrdist3

    print('\n*** Clustering %s %s chains with tcrdist3' % (len(seqs), chain_selection))

    t0 = time.time()
    tr = TCRrep(cell_df=seqs,
                organism='human',
                chains=chain,
                infer_all_genes=True,
                infer_cdrs=True,
                infer_index_cols=True,
                store_all_cdr=True,
                deduplicate=False,
                compute_distances=False)

    tr.cpus = cpus
    tr.compute_sparse_rect_distances(radius=radius, chunk_size=500)
    if chain_selection == 'alpha':
        S = tr.rw_alpha
    elif chain_selection == 'beta':
        S = tr.rw_beta
    else:
        S = tr.rw_beta

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
    mapper = cluster_TCRDist_matrix(S, seqs, method = method, cpus=cpus, hyperparam=hyper)
    t3 = time.time()

    tdist = t1-t0
    tclust = t3 - t2

    df.loc[:,'cluster']=df['bio'].map(mapper)
    df=df.rename(columns=namedict)
    # df=df[~df['cluster'].isnull()]

    return df, [tdist,tclust]

from functions.tcrdist_cpp import cluster_TCRDist_matrix_cpp

def run_tcrdistcpp(input_df, chain_selection, cpus=1, method = 'DBSCAN',radius = 50, hyper = None):

    if chain_selection =='alpha':
        tcrs = [((input_df.iloc[i]['v.alpha'],None, input_df.iloc[i]['cdr3.alpha']), (None,None,None)) for i in range(len(input_df))]
    elif chain_selection =='beta':
        tcrs = [((None,None, None), (input_df.iloc[i]['v.beta'],None, input_df.iloc[i]['cdr3.beta'])) for i in range(len(input_df))]
    else:
        tcrs = [((input_df.iloc[i]['v.alpha'],None, input_df.iloc[i]['cdr3.alpha']), (input_df.iloc[i]['v.beta'],None,input_df.iloc[i]['cdr3.beta'])) for i in range(len(input_df))]

    input_df.loc[:,'bio']=tcrs
    uniques = list(set(tcrs))
    t0 = time.time()
    D2 = calc_tcrdist_matrix_cpp(uniques,'human',
                                 single_chain=chain_selection,
                                 threshold = radius)
    # print('Retention: TCRD', np.shape(D2)[0]/len(input_df) )
    t1=time.time()
    print('Distance matrix calculated in %s seconds, clustering with %s' % (t1 - t0, method))
    t2=time.time()
    mapper = cluster_TCRDist_matrix_cpp(D2, 
                                    uniques,
                                    method,
                                    cpus=cpus,
                                    hyperparam=hyper)
    t3 = time.time()
    input_df.loc[:,'cluster']=input_df['bio'].map(mapper)
    # input_df=input_df[~input_df['cluster'].isnull()]
    tdist = t1-t0
    tclust = t3 - t2
    return input_df, [tdist,tclust]
