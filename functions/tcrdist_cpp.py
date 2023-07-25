import pandas as pd
import numpy as np
import os
from os.path import exists
from pathlib import Path
import random
# import conga: you would need to point this to wherever you downloaded the repository
# this is the top level folder for the repository (ie, it should contain scripts/ conga/ examples/)
path_to_conga = './conga'
import sys
sys.path.append(path_to_conga)
import conga
import numpy as np # not needed for tcrdist
import pandas as pd # not needed for tcrdist

from conga import util
from sklearn.cluster import DBSCAN, KMeans
import time

def distmat(tcrs, ids, distances):
    idx = pd.read_csv(ids,header=None).values.tolist()
    dists = pd.read_csv(distances,header=None).values.tolist()
    arr=np.full((len(tcrs),len(tcrs)),-1)
    idx, dists = [[[int(xi) for xi in x[0].split(' ')] for x in dx] for dx in [idx,dists]]
    for i, id in enumerate(idx):
        distlist=dists[i]
        for j,ij  in enumerate(id):
            arr[i][ij]=distlist[j]
    return arr

def calc_tcrdist_matrix_cpp(
        tcrs,
        organism,
        tmpfile_prefix = None,
        verbose=False,
        single_chain=None,
        threshold=50,
):
    if tmpfile_prefix is None:
        tmpfile_prefix = Path('./tmp_tcrdists{}'.format(random.randrange(1,10000)))
    tcrs_filename = str(tmpfile_prefix) +'_tcrs.tsv'
    df = pd.DataFrame(dict(va=[x[0][0] for x in tcrs], cdr3a=[x[0][2] for x in tcrs],
                           vb=[x[1][0] for x in tcrs], cdr3b=[x[1][2] for x in tcrs]))
    df.to_csv(tcrs_filename, sep='\t', index=False)

    if single_chain:
        name = 'find_neighbors_single_chain'
    else:
        name='find_neighbors'

    if os.name == 'posix':
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) , name)
    else:
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) , '%s.exe'%(name))

    if not exists(exe):
        print('need to compile c++ exe:', exe)
        exit(1)

    db_filename = Path.joinpath( Path(util.path_to_tcrdist_cpp_db), 'tcrdist_info_{}.txt'.format( organism ) )
    if not exists(db_filename):
        print('need to create database file:', db_filename)
        exit(1)
    if single_chain:

        cmd = '{} -c {} -t {} -f {} -d {} -o {}'.format(exe, single_chain[0].upper(),str(threshold), tcrs_filename, db_filename, tmpfile_prefix)
    else:
        cmd = '{} -f {} --only_tcrdists -d {} -o {}'.format(exe, tcrs_filename, db_filename, tmpfile_prefix)

    util.run_command(cmd, verbose=verbose)

    if single_chain:
        ids = str(tmpfile_prefix) +'_nbr{}_indices.txt'.format(threshold)
        distances = str(tmpfile_prefix) +'_nbr{}_distances.txt'.format(threshold)
        D = distmat(tcrs, ids, distances)

    else:
        tcrdist_matrix_filename = str(tmpfile_prefix) +'_tcrdists.txt'

        if not exists(tcrdist_matrix_filename):
            print('find_neighbors failed, missing', tcrdist_matrix_filename)
            exit(1)
        
        D = np.loadtxt(tcrdist_matrix_filename).astype(float)

    if single_chain:
        files = [tcrs_filename, 
                 ids, 
                 distances]
    else:
        files = [tcrs_filename, tcrdist_matrix_filename]
    for filename in files:
        os.remove(filename)

    return D

from collections import Counter
from functions.util_functions import write_lines

def cluster_TCRDist_matrix_cpp(S, seqs, method,cpus=1,hyperparam=None):
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
    assert method in methods, r'Error: %s, please choose one of the following: /n %s' % (method, methods)

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
        # inertia = kmeans.inertia_
        # write_lines('results/kmeans.csv',[[hyperparam, inertia]])
    
    # cnts = Counter(labels)
    # shortlist = [l for l in cnts.keys() if cnts[l]>=2]
    return {seq: label for seq, label in zip(seqs,labels) if label !=-1}

# def get_tcrdist_subset(input_df,chain_selection):

#     info = conga.tcrdist.all_genes.all_genes['human']
#     if chain_selection in ['alpha','beta']:
#         input_df=input_df[input_df['v.%s'%(chain_selection)].isin(info.keys())]
#     else:
#         for col in ['v.alpha','v.beta']:
#             input_df=input_df[input_df[col].isin(info.keys())]
    
#     return input_df
