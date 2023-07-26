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


def distmat(tcrs, ids, distances):
    '''Extract distance matrices from tcrdsit output file
    :param tcrs: input TCRs
    :type tcrs: Pandas DataFrame
    :param ids: TCR ID file
    :type ids: str
    :param distances: TCR distances file
    :type ids: str
    :return arr: tcrdist matrix
    :rtype arr: array

    '''
    idx = pd.read_csv(ids,header=None).values.tolist()  # Get IDs
    dists = pd.read_csv(distances,header=None).values.tolist()  # Get distances
    arr=np.full((len(tcrs),len(tcrs)),-1)   # Create input array 
    idx, dists = [[[int(xi) for xi in x[0].split(' ')] for x in dx] for dx in [idx,dists]] # Parse values
    # Map TCRs to distances
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
    '''Generate tcrdist matrix with tcrdist C++
    :param tcrs: input TCRs
    :type tcrs: list of tuples
    :param organism: organism for tcrdist reference
    :type organism: str
    :param tmpfile_prefix: output name prefix
    :type tmpfile_prefix: str
    :param verbose: set verbose
    :type verbose: bool
    :param single_chain: use single chain input
    :type single_chain: bool
    :param threshold: tcrdist radius threshold
    :type threshold: int
    :return D: distance matrix
    :type D: array'''

    # Set filename
    if tmpfile_prefix is None:
        tmpfile_prefix = Path('./tmp_tcrdists{}'.format(random.randrange(1,10000)))
    tcrs_filename = str(tmpfile_prefix) +'_tcrs.tsv'
    
    # Convert tcr input to dataframe and read to file 
    df = pd.DataFrame(dict(va=[x[0][0] for x in tcrs], cdr3a=[x[0][2] for x in tcrs],
                           vb=[x[1][0] for x in tcrs], cdr3b=[x[1][2] for x in tcrs]))
    df.to_csv(tcrs_filename, sep='\t', index=False)

    # Set tcrdist C++ implementation
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
    
    # Find database file
    db_filename = Path.joinpath( Path(util.path_to_tcrdist_cpp_db), 'tcrdist_info_{}.txt'.format( organism ) )
    if not exists(db_filename):
        print('need to create database file:', db_filename)
        exit(1)
    
    # Run tcrdist C++
    if single_chain:
        cmd = '{} -c {} -t {} -f {} -d {} -o {}'.format(exe, single_chain[0].upper(),str(threshold), tcrs_filename, db_filename, tmpfile_prefix)
    else:
        cmd = '{} -f {} --only_tcrdists -d {} -o {}'.format(exe, tcrs_filename, db_filename, tmpfile_prefix)
    util.run_command(cmd, verbose=verbose)

    # Extract distance matrix 

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

    # Clear temporary files
    if single_chain:
        files = [tcrs_filename, 
                 ids, 
                 distances]
    else:
        files = [tcrs_filename, tcrdist_matrix_filename]
    for filename in files:
        os.remove(filename)

    return D

def cluster_TCRDist_matrix_cpp(S, seqs, method,cpus=1,hyperparam=None):
    '''
    Cluster distance matrix from tcrdist
    adapted from https://github.com/svalkiers/clusTCR
    :param S: Distance matrix from tcrdist
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
    
    return {seq: label for seq, label in zip(seqs,labels) if label !=-1}
