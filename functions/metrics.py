import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_fscore_support, classification_report, adjusted_mutual_info_score, balanced_accuracy_score

def get_purity(df):
    
    df.loc[:,'epitope']=df['epitope'].replace([np.nan,'',' '],'None')
    baseline = {ep:len(df[df['epitope']==ep])/len(df) for ep in df['epitope'].unique()}
    
    mapper ={}
    purity = {}
    purity_enriched = {}
    frequency = {}
    enriched ={}
    N= {}
    Nd = 0
    clusters = df['cluster'].value_counts().index
    for c in  clusters:
        t = df[df['cluster']==c]
        N[c]=len(t)
        if len(t)>0:
            ep=t['epitope'].value_counts().index[0]
            freq = {ep:(len(t[t['epitope']==ep])/len(t)) for ep in t['epitope'].unique()}
            enrich = {k:(v/baseline[k]) for k,v in freq.items()}
            ep_enrich = max(enrich, key=enrich.get)
            mapper[c]=ep
            purity[c]=len(t[t['epitope']==ep])/len(t)
            frequency[c]=freq
            enriched[c]=ep_enrich
            purity_enriched[c]=len(t[t['epitope']==ep_enrich])/len(t)
            if len(t)<=10:
                Nd+=1
    
    return {'most_frequent':mapper,
            'most_enriched':enriched,
            'purity_frequent': purity,
            'purity_enriched':purity_enriched,
            'N': N,
            'Nd': Nd}

def get_clustermetrics(df):
    # Find clustered TCRs
    sub = df[~df['cluster'].isnull()]

    # Compute overall purity metrics
    purity_dict = get_purity(sub)

    # Compute predictive metrics
    ypred = sub['cluster'].map(purity_dict['most_frequent'])
    ytrue = sub['epitope']
    try:
        accuracy = balanced_accuracy_score(ytrue,ypred,adjusted=True)
    except ZeroDivisionError:
        print('Classification failed, zero classes')
        accuracy = 0
    precision, recall, f1score, support = precision_recall_fscore_support(ytrue,ypred,average='weighted',zero_division=0) # Total scores
    ami = adjusted_mutual_info_score(ytrue,ypred)
    rep = classification_report(ytrue, ypred, output_dict=True,zero_division=0) # Scores per epitope

    # Compute epitope-specific metrics
    counts = {k:v for k,v in sub['epitope'].value_counts().reset_index().values.tolist()}

    epmetrics = {key:{ep:rep[ep][key] for ep in counts.keys()} for key in ['precision','recall','f1-score','support']}
    maincluster = {ep: sub[sub['epitope']==ep]['cluster'].value_counts().index[0] for ep in sub['epitope'].unique()}
    mosfreq=purity_dict['most_frequent']
    clusts = {ep: '-1' if ep not in list(mosfreq.values()) else [c for c in mosfreq.keys() if mosfreq[c]==ep] for ep in counts.keys()}
    puritymap = {ep: 0 if clusts[ep]=='-1' else np.mean([len(sub[(sub['cluster']==c)&(sub['epitope']==ep)])/len(sub[sub['cluster']==c]) for c in clusts[ep]]) for ep in clusts.keys()}
    retmap = {ep: len(df[(df['epitope']==ep)&(~df['cluster'].isnull())])/len(df[df['epitope']==ep]) for ep in counts.keys()}
    consistencymap = {ep: len(sub[(sub['epitope']==ep)&(sub['cluster']==maincluster[ep])])/counts[ep] for ep in counts.keys()}
    epmetrics['consistency']=consistencymap
    epmetrics['retention']=retmap
    epmetrics['purity']=puritymap
    
    
    return {'purity': np.mean(list(purity_dict['purity_frequent'].values())),           # Purity of all clusters weighted equally (frequency)
            'purity_enriched': np.mean(list(purity_dict['purity_enriched'].values())),  # Purity of clusters (enrichment)
            'retention': len(df[~df['cluster'].isnull()])/len(df), # Proportion of clustered TCRs
            'consistency': np.mean([(consistencymap[ep]*counts[ep])/len(sub) for ep in consistencymap.keys()]), # Proportion of an epitope assigned to a given cluster
            'ami':ami,
            'accuracy':accuracy,
            'precision':precision,  # Precision over all epitopes
            'recall':recall,    # Recall over all epitopes
            'f1-score':f1score, # F1 over all epitopes
            'support':support,
            'mean_clustsize': np.mean(list(purity_dict['N'].values())),
            'small_clusters': purity_dict['Nd'],
            }, epmetrics, purity_dict


def score(df, header):
    clusterscores, epscores, clusterstats = get_clustermetrics(df)
    c_scores, e_scores = [pd.DataFrame.from_dict(d,orient='index').T for d in [clusterscores, epscores]]

    head = pd.DataFrame.from_dict(header,orient='index').T

    stats  = pd.DataFrame.from_dict(clusterstats).reset_index().rename(columns={'index':'cluster'})
    c_scores=pd.concat([head,c_scores],axis=1).reset_index()

    cols = ['epitope']+e_scores.columns.tolist()
    e_scores = e_scores.reset_index()
    e_scores.columns=cols
    for key in header.keys():
        e_scores[key]=[header[key]]*len(e_scores)
        stats[key]=[header[key]]*len(stats)

    return c_scores, e_scores, stats