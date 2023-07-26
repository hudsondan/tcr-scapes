import numpy as np
import os
import pandas as pd
path_to_conga = './conga'
import sys
sys.path.append(path_to_conga)
import conga

def strip_oovs(df, chain_selection):
    '''Remove out of vocabulary values from CDR3's
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection
    :type chain_selection: str
    :return df: cleaned data
    :type df: Pandas DataFrame'''

    # Permitted values
    vocab = ['A', 'C', 'D', 'E', 'F',
                    'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R',
                    'S', 'T', 'V', 'W', 'Y', 'X']
    
    # Retain only instances containing permitted values
    for chain in chain_selection:
        col = 'cdr3.%s' % (chain)
        oovs = [x for x in df[col].unique() if type(x) not in [float,np.float64] and len([aa for aa in x if aa not in vocab])!=0]
        df=df[~df[col].isin(oovs)].dropna(subset=col)

    return df

def get_tcrdist_subset(df, chain_selection):
    '''Retain only those TCRs whose V and J gene codes are contained
    in the IMGT reference
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection
    :type chain_selection: str
    :return df: cleaned data
    :type df: Pandas DataFrame
    '''

    info = conga.tcrdist.all_genes.all_genes['human']
    if chain_selection in ['alpha','beta']:
        df=df[df['v.%s'%(chain_selection)].isin(info.keys())]
    else:
        for col in ['v.alpha','v.beta']:
            df=df[df[col].isin(info.keys())]
    
    return df

def prepare(df, chain_selection, tcrd_sub = None):
    '''Preprocess input data
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection
    :type chain_selection: str
    :param tcrd_sub: subset to V/J gene codes in IMGT reference
    :type tcrd_sub: bool
    :return df: cleaned data
    :type df: Pandas DataFrame
    '''

    # Drop missing values
    for chain in chain_selection:
        for col in ['cdr3.%s'%(chain),
                    'v.%s'%(chain),
                    'j.%s'%(chain)]:
            if col in df.columns:
                df=df[~df[col].isnull()]
    
    # Remove instances whose V/J gene codes are missing in the IMGT reference set
    if tcrd_sub:
        print('Getting IMGT reference subset')
        info = conga.tcrdist.all_genes.all_genes['human']
        for chain in chain_selection:
            for gene in ['v', 'j']:
                col = '%s.%s' % (gene,chain)
                if col in df.columns:
                    df[col]=df[col].map({x:x if x in info.keys() else x.split('*')[0]+'*01' if x.split('*')[0]+'*01' in info.keys() else np.nan for x in df[col].unique()})
                    df=df[~df[col].isnull()]

    return df

def preprocess(data, chain_selection, n=0, paired=None, tcrd_sub=None):
    '''Process input data for clustering
    :param data: input data
    :type data: Pandas DataFrame
    :param chain_selection: TCR chain selection
    :type chain_selection: str
    :param n: minimum number of TCRs per epitope
    :type n: int
    :param paired: subset data for paired instances only
    :type paired: bool
    :param tcrd_sub: subset to V/J gene codes in IMGT reference
    :type tcrd_sub: bool
    :return df: cleaned data
    :type df: Pandas DataFrame'''

    # Retain only paired chain instances
    if paired:
        print('Retaining paired chain instances')
        data = data[data['pairing']=='paired']
    # Retain only instances with epitope linkage
    df=data[data['epitope']!='None']
    
    # Retain instances with length > 5 (tcrdist c++ requirement)
    df=df[(df['length_alpha']>5)&(df['length_beta']>5)]    
    df = strip_oovs(df,chain_selection) # Remove out of vocab values

    cols = [x for x in ['cdr3.alpha',
                        'v.alpha',
                        'j.alpha',
                        'v.alpha_clean_level_3',
                        'j.alpha_clean_level_3',
                        'cdr3.beta',
                        'v.beta',
                        'j.beta',
                        'v.beta_clean_level_3',
                        'j.beta_clean_level_3',
                        'subject',
                        'epitope',
                        'pairing',
                        'dataset'] if x in df.columns]
    df=df[cols]
    for c in ['v.alpha',
              'j.alpha',
              'v.beta',
              'j.beta']:
        if (c in cols) & (c+'_clean_level_3' in cols):
            df=df.drop(labels=c,axis='columns')
            df=df.rename(columns={c+'_clean_level_3':c})

    # Filter for missing values and (optionally) IMGT gene codes
    df = prepare(df,
                 chain_selection, 
                 tcrd_sub=tcrd_sub)
    
    # Retain well-represented epitopes
    if n>0:
        print('Retaining epitopes with ≥ %s representatives'%(n))
        counts= df['epitope'].value_counts().reset_index()
        eps=counts[counts['count']>=int(n)]['epitope'].values.tolist()
        df=df[df['epitope'].isin(eps)]
    return df

def get_bio(df, chain_selection, add_j=True):
    '''Get bioidentity for TCR instances
    :param df: input data
    :type df: Pandas DataFrame
    :param chain_selection: TCR chain selection
    :type chain_selection: str
    :param add_j: include J gene codes
    :type add_j: bool
    :return df: updated dataframe
    :rtype df: Pandas DataFrame'''
    
    if add_j:
        print('Getting Bio-ID on V, CDR3, J')
        cols = ['%s.%s'%(x, chain_selection) for x in ['v','cdr3','j']]

    else:
        print('Getting Bio-ID on V, CDR3')
        cols = ['%s.%s'%(x, chain_selection) for x in ['v','cdr3']]

    df.loc[:,'bio']=['-'.join(x) for x in df[cols].values.tolist()]
    
    return df

def get_olga(chain_selection, n):
    '''Generate random TCRs with OLGA (Sethna et al 2019)
    :param chain_selection: TCR chain
    :type chain_selection: str
    :param n: # of sequences to generate
    :type n: int
    :return tr: TCRs
    :rtype tr: Pandas DataFrame'''
    
    root='functions/olga/'
    path = os.path.join(root,'olga_TR_{}.csv'.format(chain_selection))
    os.system('python functions/olga/generate_sequences.py --human_T_{} -n {} -o {}'.format(chain_selection, n,path))
    tr =pd.read_csv(path,header=None)
    tr.columns=['cdr3.%s.nt'%(chain_selection),'cdr3.%s'%(chain_selection),'v.%s'%(chain_selection),'j.%s'%(chain_selection)]
    tr['epitope']=['DECOY']*len(tr)
    os.system('rm {}'.format(path))
    return tr