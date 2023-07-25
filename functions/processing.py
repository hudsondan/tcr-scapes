from functions.util_functions import load_pickle
import numpy as np
import os
import pandas as pd
path_to_conga = './conga'
import sys
sys.path.append(path_to_conga)
import conga

def strip_oovs(df, chain_selection):
    vocab = ['A', 'C', 'D', 'E', 'F',
                    'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R',
                    'S', 'T', 'V', 'W', 'Y', 'X']
    for chain in chain_selection:
        col = 'cdr3.%s' % (chain)
        oovs = [x for x in df[col].unique() if type(x) not in [float,np.float64] and len([aa for aa in x if aa not in vocab])!=0]
        df=df[~df[col].isin(oovs)].dropna(subset=col)
    return df

def drop_multis(df,pairing):
    if pairing in ['alpha','paired']:
        multia = [tcr for tcr in df['cdr3.alpha'].unique() if len(df[df['cdr3.alpha']==tcr]['epitope'].unique())>1]
        df = df[~df['cdr3.alpha'].isin(multia)]

    if pairing in ['beta','paired']:
        multib = [tcr for tcr in df['cdr3.beta'].unique() if len(df[df['cdr3.beta']==tcr]['epitope'].unique())>1]
        df = df[~df['cdr3.beta'].isin(multib)]
    return df

def get_tcrdist_subset(input_df,chain_selection):

    info = conga.tcrdist.all_genes.all_genes['human']
    if chain_selection in ['alpha','beta']:
        input_df=input_df[input_df['v.%s'%(chain_selection)].isin(info.keys())]
    else:
        for col in ['v.alpha','v.beta']:
            input_df=input_df[input_df[col].isin(info.keys())]
    
    return input_df

def prepare(df, chain_selection, tcrd_sub = None):
    for chain in chain_selection:
        for col in ['cdr3.%s'%(chain),
                    'v.%s'%(chain),
                    'j.%s'%(chain)]:
            if col in df.columns:
                df=df[~df[col].isnull()]
    if tcrd_sub:
        print('Getting tcrdist subset')
        info = conga.tcrdist.all_genes.all_genes['human']
        for chain in chain_selection:
            for gene in ['v', 'j']:
                col = '%s.%s' % (gene,chain)
                if col in df.columns:
                    df[col]=df[col].map({x:x if x in info.keys() else x.split('*')[0]+'*01' if x.split('*')[0]+'*01' in info.keys() else np.nan for x in df[col].unique()})
                    df=df[~df[col].isnull()]

    return df

def generate_data(data, chain_selection, n=0, paired=None, tcrd_sub=None):

    if paired:
        print('Retaining paired chain instances')
        data = data[data['pairing']=='paired']
    df=data[data['epitope']!='None']
    df=df[(df['length_alpha']>5)&(df['length_beta']>5)]    
    df = strip_oovs(df,chain_selection)
    cols = ['cdr3.alpha',
            'v.alpha_clean_level_3',
            'j.alpha_clean_level_3',
            'cdr3.beta',
            'v.beta_clean_level_3',
            'j.beta_clean_level_3',
            'subject',
            'epitope',
            'pairing',
            'dataset']
    df=df[cols]
    df = df.rename(columns={'v.alpha_clean_level_3':'v.alpha',
                            'j.alpha_clean_level_3':'j.alpha',
                            'v.beta_clean_level_3':'v.beta',
                            'j.beta_clean_level_3':'j.beta',
                            })
    df = prepare(df,
                 chain_selection, 
                 tcrd_sub=tcrd_sub)
    if n>0:
        print('Retaining epitopes with ≥ %s representatives'%(n))
        counts= df['epitope'].value_counts().reset_index()
        eps=counts[counts['count']>=int(n)]['epitope'].values.tolist()
        df=df[df['epitope'].isin(eps)]
    return df

def get_bio(df,pairing,add_j=True):
    
    if add_j:
        print('Getting Bio-ID on V, CDR3, J')
        cols = ['%s.%s'%(x, pairing) for x in ['v','cdr3','j']]

    else:
        print('Getting Bio-ID on V, CDR3')
        cols = ['%s.%s'%(x, pairing) for x in ['v','cdr3']]

    df.loc[:,'bio']=['-'.join(x) for x in df[cols].values.tolist()]
    
    return df

def get_olga(chain_selection, n):
    root='functions/olga/'
    path = os.path.join(root,'olga_TR_{}.csv'.format(chain_selection))
    os.system('python functions/olga/generate_sequences.py --human_T_{} -n {} -o {}'.format(chain_selection, n,path))
    tr =pd.read_csv(path,header=None)
    tr.columns=['cdr3.%s.nt'%(chain_selection),'cdr3.%s'%(chain_selection),'v.%s'%(chain_selection),'j.%s'%(chain_selection)]
    tr['epitope']=['DECOY']*len(tr)
    os.system('rm {}'.format(path))
    return tr