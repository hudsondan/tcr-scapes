import numpy as np
import pandas as pd

def get_pairing(df):
    if 'pairing' not in df.columns:
        p = [[False if (type(x)==float)|(type(x)==np.float64) else len(x)!=0 for x in df['cdr3.%s'%(chain)].values] for chain in ['alpha','beta']]
        pairing = [(p[0][i],p[1][i]) for i in range(len(df))]
        df['pairing']=['paired' if x ==(True,True) else 'alpha' if x==(True,False) else 'beta' for x in pairing]
    return df

def clip(seq, pad = 50):
    '''Zero pad TCR sequences to a common length
    :param seq: input sequence
    :type seq: str
    :param pad: pad length
    :type pad: int
    :return: zero-padded sequence
    :rtype: str
    '''
    if type(seq) not in [float,np.float64]:
        seq=seq.rstrip()
    padlen = pad if type(seq) in [float,np.float64] else 0 if len(seq)>pad else pad-len(seq)
    return'X'*pad if type(seq)in [float,np.float64] else seq[:pad]+'X'*padlen

def clean_vdj_input(df):
    for col in ['v.alpha','j.alpha','v.beta','j.beta']:
        if col in df.columns:
            df[col+'_clean_level_1']=df[col].apply(lambda x: 'None' if type(x) in [float,np.float64] else ''.join(x.split('*')[0].split('-')[0].split(' ')))
            df[col+'_clean_level_2']=df[col].apply(lambda x: 'None' if type(x) in [float,np.float64] else ''.join(x.split('*')[0]))
            df[col+'_clean_level_3']=df[col].apply(lambda x: 'None' if type(x) in [float,np.float64] else ''.join(x.split(':')[0]))
    return df

def get_unique(df):
    '''Generate a single sequence representative of cdr3 alpha and beta
    :param df: input dataframe
    :type df: DataFrame
    :return df: output dataframe
    :rtype df: DataFrame'''
    
    for col in ['cdr3.alpha','cdr3.beta']:
        df[col+'_processed'] = [clip(x) for x in df[col].values]
    
    if (('v.alpha_clean_level_1' not in df.columns)|('v.beta_clean_level_1' not in df.columns)):
        df = clean_vdj(df)
    df['bio_identity_alpha']=df['cdr3.alpha_processed']+'_'+df['v.alpha_clean_level_1']+'_'+df['j.alpha_clean_level_1']
    df['bio_identity_beta']=df['cdr3.beta_processed']+'_'+df['v.beta_clean_level_1']+'_'+df['j.beta_clean_level_1']
    df['unique_cdr3']=df['bio_identity_alpha']+'_'+df['bio_identity_beta']
    if 'antigen.epitope_clean' in df.columns:
        df['unique_cdr3_ep']=df['bio_identity_alpha']+'_'+df['bio_identity_beta']+'_'+df['antigen.epitope_clean']

    return df

def clean_species(df):
    if 'antigen.species' in df.columns:
        df['antigen.species_clean']=df['antigen.species'].replace({
            np.nan: 'None',
            'Acidiphilium cryptum (Acidiphilum cryptum)':'Acidiphilium cryptum',
            'AKR (endogenous) murine leukemia virus (AKR murine leukemia virus)': 'AKR-MLV',
            'Alkalihalobacillus clausii KSM-K16':'Alkalihalobacillus clausii',
            'Arabidopsis thaliana (mouse-ear cress)':'Arabidopsis thaliana',
            'Arachis hypogaea (peanut)':'Arachis hypogaea',
            'Bacillus [genus]': 'Bacillus',
            'Bat mastadenovirus G':'Bat mastadenovirus',
            'Blastococcus saxobsidens DD2 (Blastococcus saxobsidens str. DD2)':'Blastococcus saxobsidens',
            'Borreliella burgdorferi (Lyme disease spirochete)':'Borreliella burgdorferi',
            'Bos taurus (bovine)': 'Bos taurus',
            'Chlorobium chlorochromatii (epibiont of the phototrophic consortium "Chlorochromatium aggregatum")': 'Chlorobium chlorochromatii',
            'Coxsackievirus B4 (Coxsackie B virus type 4)':'Coxsackievirus',
            'dengue virus type I (Dengue virus 1)':'DENV',
            'Dengue virus type 2 (Dengue virus 2)':'DENV',
            'Dengue virus type 3 (Dengue virus serotype 3)':'DENV',
            'DENV-1':'DENV',
            'DENV-2':'DENV',
            'DENV-3/4':'DENV',
            'DENV1':'DENV',
            'DENV2':'DENV',
            'DENV3/4':'DENV',
            'Dengue virus 1 (Dengue virus type 1)': 'DENV',
            'Dengue virus 2 (Dengue virus type 2)': 'DENV',
            'Dengue virus 3 (Dengue virus serotype 3)':'DENV',
            'Desulfotignum phosphitoxidans DSM 13687':'Desulfotignum phosphitoxidans',
            'Encephalitozoon romaleae SJ-2008':'Encephalitozoon romaleae',
            'E.Coli': 'E. coli',
            'H1N1 subtype (H1N1)':'H1N1',
            'HSV2':'HSV-2',
            'HTLV1':'HTLV-1',
            'HTLV-1-chronic':'HTLV-1',
            'Hepatitis B virus (Human hepatitis B virus)':'HBV',
            'Hepatitis C virus':'HCV',
            'Hepatitis C virus subtype 3a (Hepatitis C virus 3a)': 'HCV',
            'Hepacivirus C':'HCV',
            'Hepatitis C virus (isolate BK)':'HCV',
            'Hepatitis E virus type 3 (Hepatitis E virus genotype 3)':'HEV',
            'HIV-1':'HIV',
            'HIV-1 M:B_HXB2R (Human immunodeficiency virus type 1 (HXB2 ISOLATE))':'HIV',
            'Human herpesvirus 3 (Varicella-zoster virus)':'HHV',
            'Human herpesvirus 5 strain AD169 (Human cytomegalovirus (strain AD169))':'HHV',
            'Human coronavirus OC43 (Human coronavirus (strain OC43))':'SARS-CoV-1',
            'Human coronavirus NL63 (Coronavirus NL63)':'SARS-CoV-1',
            'Human herpesvirus 2':'HHV',
            'Homo sapiens (human)': 'Homo sapiens',
            'HomoSapiens': 'Homo sapiens',
            'Hordeum vulgare (barley)': 'Hordeum vulgare',
            'Human T-cell leukemia virus type I (Human T cell leukemia virus type 1)':'HTLV-1',
            'Human herpesvirus 1': 'HHV-1',
            'Human herpesvirus 4 (Epstein Barr virus)':'EBV',
            'Human herpesvirus 5 (Human cytomegalovirus)': 'CMV',
            'Human immunodeficiency virus 1 (human immunodeficiency virus 1 HIV-1)':'HIV',
            'Human papillomavirus type 16 (Human papilloma virus type 16)':'HPV',
            'InfluenzaA': 'Influenza',
            'Influenza A virus':'Influenza',
            'Influenza A virus (A/Canterbury/236/2005(H3N2))':'Influenza',
            'Influenza A virus (A/Memphis/4/1973(H3N2)) (Influenza A virus (A/Memphis/4/73(H3N2)))':'Influenza',
            'Influenza A virus (A/bar-headed goose/Qinghai/3/2005(H5N1)) (Influenza A virus (A/bar headed goose/Qinghai/3/2005(H5N1)))':'Influenza',
            'Influenza A virus (A/chicken/Anhui/1/1998(H9N2)) (Influenza A virus (A/Chicken/Anhui/1/98(H9N2)))':'Influenza',
            'Influenza A virus (A/California/07/2009(H1N1))':'Influenza',
            'Influenza A virus (A/Puerto Rico/8/1934(H1N1)) (Influenza A virus (A/PR 8/34 (H1N1)))':'Influenza',
            'Influenza A virus (A/California/04/2009(H1N1))':'Influenza',
            'JC polyomavirus (Human polyomavirus (type JC))': 'JC virus',
            'Kitasatospora setae KM-6054 (Kitasatospora setae str. KM-6054)':'Kitasatospora setae',
            'Legionella longbeachae NSW150 (Legionella longbeachae str. NSW150)':'Legionella longbeachae',
            'M.tuberculosis': 'Mycobacterium tuberculosis',
            'M. tuberculosis': 'Mycobacterium tuberculosis',
            'Macrophomina phaseolina (charcoal rot)': 'Macrophomina phaseolina',
            'Merkel cell polyomavirus (MCPyV isolate R17b)':'MCPyV',
            'Mtb': 'Mycobacterium tuberculosis',
            'Mus musculus (mouse)': 'Mus musculus',
            'Mycobacterium tuberculosis H37Rv (Mycobacterium tuberculosis str. H37Rv)': 'Mycobacterium tuberculosis',
            'Oryctolagus cuniculus (rabbit)':'Oryctolagus cuniculus',
            'Oryzias latipes (Japanese medaka)': 'Oryzias latipes',
            'PseudomonasFluorescens':  'Pseudomonas fluorescens',
            'PseudomonasAeruginosa': 'Pseudomonas aeruginosa',
            'Rhodococcus sp. AW25M09':'Rhodococcus sp.',
            'SARS coronavirus BJ01':'SARS-CoV-1',
            'SARS coronavirus Tor2 (Severe acute respiratory syndrome-related coronavirus Tor2)':'SARS-CoV-1',
            'SARS-CoV1':'SARS-CoV-1',
            'SARS-CoV2':'SARS-CoV-2',
            "Saccharomyces cerevisiae (baker's yeast)": 'Saccharomyces cerevisiae',
            'SaccharomycesCerevisiae':'Saccharomyces cerevisiae',
            'Salmonella enterica subsp. enterica serovar Typhimurium':'Salmonella enterica',
            'Schinkia azotoformans (Bacillus azotiformans)':'Schinkia azotoformans',
            'Shigella dysenteriae serotype 1':'Shigella dysenteriae',
            'StreptomycesKanamyceticus': 'Streptomyces kanamyceticus',
            'SelaginellaMoellendorffii':'Selaginella moellendorffii',
            'Sulfurovum sp. NBC37-1':'Sulfurovum sp.',
            'Trichosporon asahii var. asahii CBS 2479': 'Trichosporon asahii',
            'Triticum aestivum (Canadian hard winter wheat)':'Triticum aestivum',
            'TriticumAestivum':'Triticum aestivum',
            'Yellow fever virus (Flavivirus febricis)': 'YFV',
            'YellowFever':'YFV',
            'Volvox carteri f. nagariensis':'Volvox carteri',

        })
    else:
        df['antigen.species_clean']=['None']*len(df)
    return df

def clean_gene(df):
    if 'antigen.gene' in df.columns:
        df['antigen.gene_clean']=df['antigen.gene'].replace({
            
            '55 kDa immediate-early protein 1':'EI1',
            '45 kDa immediate-early protein 2':'IE2',
            '70KFP (U1Ð70 kDa fusionprotein)':'70KFP',
            'BMLF1 protein':'BMLF1',
            'C.albicans':'C. albicans',
            'Ca2-indepen-Plip-A2':'iPLA2',
            'calcium-independent phospholipase A2 [Homo sapiens]':'iPLA2',
            'DQ2.5-glia-?1':'DQ2.5-glia',
            'DQ2.5-glia-?1a':'DQ2.5-glia',
            'DQ2.5-glia-?1a|DQ2.5-glia-?1':'DQ2.5-glia',
            'DQ2.5-glia-?1|DQ2.5-glia-?2':'DQ2.5-glia',
            'DQ2.5-glia-?2':'DQ2.5-glia',
            'DQ2.5-glia-?2|DQ2.5-glia-?2':'DQ2.5-glia',
            'DQ2.5-glia-?2|DQ2.5-glia-?4c|DQ2.5-glia-?4d|DQ2.5-glia-?2':'DQ2.5-glia',
            'DQ2.5-glia-?4a|DQ2.5-glia-?4b|DQ2.5-glia-?1':'DQ2.5-glia',
            'DQ2.5-glia-?4c|DQ2.5-glia-?2':'DQ2.5-glia',
            'gag polyprotein':'Gag',
            'gag protein':'Gag',
            'GagpolyproteinRQ13':'Gag',
            'EBNA-3A nuclear protein':'EBNA-3A',
            'EBNA3A nuclear protein':'EBNA-3A',
            'EBNA-3B nuclear protein':'EBNA-3B',
            'EBNA-3C latent protein':'EBNA-3C',
            'EBNA3C latent protein':'EBNA-3C',
            'envelope':'Env',
            'envelope,ORF1ab':'Env/ORF1ab',
            'hepatitis B surface antigen':'HBSag',
            'human T-cell leukemia virus(HTLV-1)':'HTLV-1',
            'HTLV-1 Tax11-19':'HTLV-1',
            'Human T-cell lymphotropic virus type 1 (HTLV-1) ':'HTLV-1',
            'insulin precursor A1-15':'Insulin',
            'Latent membrane protein 2':'LMP2',
            'Low molecular weight glutenin subunit precursor':'LMW-GS',
            'Matrixprotein(M1)': 'M-protein',
            'Matrix protein 1':'M-protein',
            'Melan-A/MART-1':'MART-1',
            'Melan-AA27L':'MART-1',
            'MelanA/MART1':'MART-1',
            'membrane glycoprotein':'M-gp',
            'membrane glycoprotein,surface glycoprotein':'M-gp/S-gp',
            'myelin basic protein':'MBP',
            'M.tuberculosis':'Mtb',
            'Nef138-10(2F)':'Nef',
            'Nef138-10(wt)':'Nef',
            'Nef138-11(wt)':'Nef',
            'Nef138-12(wt)':'Nef',
            'Nef138-13(wt)':'Nef',
            'Nef138-14(wt)':'Nef',
            'Nef138-15(wt)':'Nef',
            'Nef138-16(wt)':'Nef',
            'non-structural protein NS4b':'NS4b',
            'nuclear antigen EBNA-3':'EBNA-3',
            'nuclear antigen EBNA1':'EBNA-1',
            'nucleocapsid phosphoprotein':'N-pp',
            '-':'None',
            'ORF1ab,ORF3a':'ORF1ab/ORF3a',
            'ORF1ab,surface glycoprotein':'ORF1ab/S-gp',
            'porin OmpC [Escherichia coli]':'OmpC',
            'Protein Nef':'Nef',
            'recombinant hepatitis B surface antigen':'HBSag',
            'regulatory protein IE1': 'IE1',
            'regulatory protein IE1 [Human betaherpesvirus 5]': 'IE1',
            'surface glycoprotein':'S-gp',
            'Tax-1':'Tax',
            'TetanusToxoid':'TT',
            'transcriptional activator Tax':'Tax',
            'zinc transporter 8': 'ZNT8',
            'zinc transporter 8 isoform a':'ZNT8',
            'nan':'None',
            
    })
    else:
        df['antigen.gene_clean']=['None']*len(df)
    return df

def clean_vdj(df):
    for col in ['v.alpha','j.alpha','v.beta','d.beta','j.beta']:
        if col in df.columns:
            mapper = pd.read_csv('Data/vdj_assignments/%s.csv'%(col))
            mapper1 = {mapper.iloc[i][col]:'None' if type(mapper.iloc[i]['%s_clean'%(col)]) in [float,np.float64] else mapper.iloc[i]['%s_clean'%(col)].split('*')[0].split('-')[0] for i in range(len(mapper))}
            mapper2 = {mapper.iloc[i][col]:'None' if type(mapper.iloc[i]['%s_clean'%(col)]) in [float,np.float64] else mapper.iloc[i]['%s_clean'%(col)].split('*')[0] for i in range(len(mapper))}
            mapper3 = {mapper.iloc[i][col]:'None' if type(mapper.iloc[i]['%s_clean'%(col)]) in [float,np.float64] else mapper.iloc[i]['%s_clean'%(col)] for i in range(len(mapper))}
            df['%s_clean_level_1'%(col)]=df[col].replace(mapper1).replace(['', ' ',np.nan],'None')
            df['%s_clean_level_2'%(col)]=df[col].replace(mapper2).replace(['', ' ',np.nan],'None')
            df['%s_clean_level_3'%(col)]=df[col].replace(mapper3).replace(['', ' ',np.nan],'None')

    return df

def get_mapper(mapper):
    mapper['input_mapped_level_1']=mapper['input_mapped'].apply(lambda x: 'None' if type(x)==float else x.strip('w').split('*')[0])
    mapper['input_mapped_level_2']=mapper['input_mapped'].apply(lambda x: 'None' if type(x)==float else x.strip('w').split(':')[0])
    mapper['input_mapped_level_3']=mapper['input_mapped'].apply(lambda x: 'None' if type(x)==float else ':'.join(x.strip('w').split(':')[:2]))
    mapper['input_mapped_2_level_1']=mapper['input_mapped_2'].apply(lambda x: 'None' if type(x)==float else x.strip('w').split('*')[0])
    mapper['input_mapped_2_level_2']=mapper['input_mapped_2'].apply(lambda x: 'None' if type(x)==float else x.strip('w').split(':')[0])
    mapper['input_mapped_2_level_3']=mapper['input_mapped_2'].apply(lambda x: 'None' if type(x)==float else ':'.join(x.strip('w').split(':')[:2]))
    mapper_dict = {}
    for i in range(len(mapper)):
        if type(mapper.iloc[i]['input_mapped_2'])==float:
            mapper_dict[mapper.iloc[i]['input']]=mapper.iloc[i]['input_mapped_level_2']
        else:
            mapper_dict[mapper.iloc[i]['input']]=mapper.iloc[i]['input_mapped_2_level_2']
    return mapper_dict
    
def get_subset(x,s):
    if type(x)==float:
        return 'None'
    else:
        l = list(set([xi for xi in x if s in xi]))
        if len(l)==0:
            return 'None'
        else:
            return '-'.join(l)

def get_hlas():
    return sorted([
                   'hla.type',
                   'HLA-A',
                    'HLA-A_1',
                    'HLA-B',
                    'HLA-B_1',
                    'HLA-C',
                    'HLA-C_1',
                    'DPA1',
                    'DPA1_1',
                    'DPB1',
                    'DPB1_1',
                    'DQA1',
                    'DQA1_1',
                    'DQB1',
                    'DQB1_1',
                    'DRB1',
                    'DRB1_1',
                    'DRB3',
                    'DRB3_1',
                    'DRB4',
                    'DRB4_1',
                    'DRB5',
                    'DRB5_1',
                    'HLA-E',
                    'CD1',
                    'MR1',
                    'HLA-E_1', 
                    ]), ['HLA-A',
                    'HLA-B',
                    'HLA-C',
                    'HLA-DPA1',
                    'HLA-DPB1',
                    'HLA-DQA1',
                    'HLA-DQB1',
                    'HLA-DRA1',
                    'HLA-DRB1',
                    'HLA-DRB3',
                    'HLA-DRB4',
                    'HLA-DRB5',
                    'HLA-E',
                    'CD1',
                    'MR1',
                    ]

def hla_map(df):
    hla_cols, shortlist = get_hlas()
    mapper_dict = get_mapper(pd.read_csv('Data/HLA_mapper_backup2.csv',index_col=0,low_memory=False))
    if 'mhc' in df.columns:
        df['mhc_mapped']=df['mhc'].replace(mapper_dict).replace(np.nan,'None')
    null = ['None',np.nan, '']
    out = []
    hla_cols = [c for c in hla_cols if c in df.columns]
    try:
        for ref_inp in df[hla_cols].values.tolist():
            ref_inp = [x for x in ref_inp if x not in null and type(x) !=float]
            if len(ref_inp)!=0:

                hlas = np.concatenate([x.split(',') for x in list(set(ref_inp)) if x not in ['None',np.nan, '']])
                hlas = np.concatenate([x.split(';') for x in hlas]).tolist()
                hlas = np.concatenate([x.split('/') for x in hlas]).tolist()
                try:
                    out.append(sorted([mapper_dict[x] for x in hlas if x !='']))
                except KeyError:
                    out2 = []
                    for x in hlas:
                        if x!='':
                            try:
                                try:
                                    out2.append(mapper_dict[x])
                                except KeyError:
                                    out2.append(mapper_dict[x.lstrip()])
                            except KeyError:
                                continue
                                
                    out.append(sorted(out2))
            else:
                out.append('None')
        df['hla.type_mapped']=out
        for s in shortlist:
            df[s+'_mapped']=df['hla.type_mapped'].apply(lambda x:get_subset(x,s))
    except KeyError:
        print('HLA mapping skipped: no corresponding columns')
        
    return df

def label_epitopes(df):
    '''Combine species, gene and epitope information
    :param df: input dataframe
    :type df: DataFrame
    :return df: output dataframe
    :rtype df: DataFrame
    '''
    if (('antigen.species_clean' in df.columns)&('antigen.gene_clean' in df.columns)&('antigen.epitope_clean' in df.columns)):
        df['species_gene_epitope']=['_'.join(x) for x in df[['antigen.species_clean','antigen.gene_clean','antigen.epitope_clean']].values.tolist()]
    return df


def replace_str(seq):
    vocab = ['A', 'C', 'D', 'E', 'F',
                'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'P', 'Q', 'R',
                'S', 'T', 'V', 'W', 'Y', 'X']
    missing = list(set([x for x in seq if x not in vocab]))
    for char in missing:
        seq = seq.replace(char,'X')
    return seq.upper()

def get_length(df):

   for chain in ['alpha','beta']:
      if 'length_%s'%(chain) not in df.columns:
        df['length_%s'%(chain)]=df['cdr3.%s'%(chain)].apply(lambda x:0 if ((type(x)==float)|(x=='None')) else len(x.strip('X')))
   df['length'] = np.int64([np.max(x) if ((x[0]==0)|(x[1]==0))  else np.mean(x) for x in df[['length_alpha','length_beta']].values.tolist()])
   return df

def clean(df):
    for col in ['cdr3.alpha','cdr3.beta']:
        df[col]=[x if type(x)  in [np.float64,float] else x.rstrip() for x in df[col].values]
    for col in ['v.alpha','j.alpha','v.beta','d.beta','j.beta']:
        if col not in df.columns:
            df[col]=['None']*len(df)
    df=get_length(df)

    for col in ['antigen.species','antigen.gene','antigen.epitope']:
        if col in df.columns:
            df[col]=df[col].replace({'':'None',
                                    ' ':'None',
                                    '-':'None',
                                    np.nan:'None',
                                    'NA':'None'})
        else:
            df[col]=['None']*len(df)
    

    mapper = {x: x.split('+')[0].rstrip() for x in [x for x in df['antigen.epitope'].unique() if '+' in x]}
    df['antigen.epitope_clean']=df['antigen.epitope'].replace(mapper)


    return hla_map(get_unique(label_epitopes(clean_gene(clean_species(clean_vdj(get_pairing(df)))))))

def split_paired(df):
    '''Split paired sequences into their constituent alpha and beta chains
    :param df: input dataframe
    :type df: DataFrame
    :return: output dataframe
    :rtype: DataFrame'''
    
    sub = df[df['pairing']=='paired']
    
    alpha=sub.copy()
    alpha['cdr3.beta']=[np.nan] * len(alpha)
    alpha['v.beta']=[np.nan] * len(alpha)
    alpha['d.beta']=[np.nan] * len(alpha)
    alpha['j.beta']=[np.nan] * len(alpha)
    alpha['pairing']=['alpha']*len(alpha)

    beta=sub.copy()
    beta['cdr3.alpha']=[np.nan] * len(beta)
    beta['v.alpha']=[np.nan] * len(beta)
    beta['j.alpha']=[np.nan] * len(beta)
    beta['pairing']=['beta']*len(beta)

    return pd.concat([df,alpha,beta])

def read_process(df_file,cols=None):
    '''Load and preprocess input dataframes
    :param df: input dataframe
    :type df: DataFrame
    :param cols: columns names:
    :type cols: list
    :return df: output dataframe
    :rtype df: DataFrame'''
    
    df = pd.read_csv(df_file,index_col=0,low_memory=False)
    df=df.reset_index(drop=False)
    if not cols:
        cols=df.columns
    
    for col in ['cdr3.alpha','cdr3.beta']:
        if col not in cols:
            print('Adding ',col)
            df[col]=[np.nan]*len(df)

    # Get chain pairing
    if 'pairing' not in df.columns:
        df=get_pairing(df)

    # Replace missing values
    for col in cols:
        if col not in ['cdr3.alpha','cdr3.beta']:
            if col in df.columns:
                df[col]=df[col].replace({'-':'None',
                            np.nan: 'None',})
            else:
                df[col]=['None']*len(df)
    
    if (('v.alpha_clean_level_1' not in df.columns)|('j.alpha_clean_level_1' not in df.columns)|('v.beta_clean_level_1' not in df.columns)|('j.beta_clean_level_1' not in df.columns)):
        df = clean_vdj_input(df)
    
    df=clean(df)

    return df