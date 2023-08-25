import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_fscore_support, classification_report, adjusted_mutual_info_score, balanced_accuracy_score
from scipy.stats import hmean
from sklearn.metrics import multilabel_confusion_matrix


def get_purity(df):
    '''Compute cluster purity
    :param df: input data
    :type df: pandas DataFrame
    :return: total metrics, epitope-secific metrics and cluster statistics
    :rtype: dict'''

    df.loc[:,'epitope']=df['epitope'].replace([np.nan,'',' '],'None')   # Recode NaNs
    baseline = {ep:len(df[df['epitope']==ep])/len(df) for ep in df['epitope'].unique()} # Get baseline frequence for eadch epitope
    
    mapper ={}
    purity = {}
    purity_enriched = {}
    frequency = {}
    enriched ={}
    N= {}
    Nd = 0
    clusters = df['cluster'].value_counts().index   # Rank clusters
    for c in  clusters:
        t = df[df['cluster']==c]    
        N[c]=len(t)
        if len(t)>0:
            ep=t['epitope'].value_counts().index[0] # Most frequent epitope
            freq = {ep:(len(t[t['epitope']==ep])/len(t)) for ep in t['epitope'].unique()}   # Get frequency of each epitope 
            enrich = {k:(v/baseline[k]) for k,v in freq.items()}    # Compute enrichment of each epitope vs. baseline
            ep_enrich = max(enrich, key=enrich.get) # Find the most enriched
            mapper[c]=ep    # Record the most frequent epitope
            purity[c]=len(t[t['epitope']==ep])/len(t)   # Compute purity for the most frequent epitope
            frequency[c]=freq   # Record the most frequent epitope
            enriched[c]=ep_enrich # Record the most enriched epitope 
            purity_enriched[c]=len(t[t['epitope']==ep_enrich])/len(t)   # Compute purity for the most enriched epitope
            if len(t)<=10:
                Nd+=1   # Find the number of clusters with 10 or fewer members
    
    return {'most_frequent':mapper,
            'most_enriched':enriched,
            'purity_frequent': purity,
            'purity_enriched':purity_enriched,
            'N': N,
            'Nd': Nd}

def precision_recall_fscore_include_uncalled(df, ytrue, ypred, dryrun=False, debug=False):
    """
    Just adding None to ypred and a dummy value to ytrue is not supported
    We could derive fp, tp etc from ypred and ytrue manually through comparison,
    but below we derive it from multilabel_confusion_matrix and also need
    to implement average='weighted'
    """

    if dryrun:
        print("DRYRUN: not actually adding missing FNs. Values should be identical with original. Only useful for debugging")

    labels = np.unique(ytrue)
    precisions = []# per epitope
    recalls = []# per epitope
    weights = []# per epitope
    fn_per_epi = df[df['cluster'].isnull()]['epitope'].value_counts()# previously ignored FNs
    for (i, cm) in enumerate(multilabel_confusion_matrix(ytrue, ypred, labels=labels)):# per epitope
        # for order see https://scikit-learn.org/stable/modules/generated/sklearn.metrics.multilabel_confusion_matrix.html
        tn = cm[0][0]
        fn = cm[1][0]
        tp = cm[1][1]
        fp = cm[0][1]
        lbl = labels[i]

        if not dryrun:
            missing_fn = fn_per_epi.get(lbl, 0)# for random there likely won't be any missing
            fn += missing_fn
            if debug:
                print(f"DEBUG {lbl} Adding {missing_fn} to fn")
                print(f"DEBUG {lbl} tp={tp} fp={fp} fn={fn} tn={tn}")

        if tp+fp == 0:
            precision = 0.0
        else:
            precision = tp/(tp+fp)

        if tp+fn == 0:
            recall = 0.0
        else:
            recall = tp/(tp+fn)

        # weighting='average'
        support = sum(ytrue == lbl)
        w = support/ytrue.shape[0]
        weights.append(w)
        precision *= w
        recall *= w
        #if debug:
        #    print(f"{lbl} after weighting precision={precision} recall={recall}")

        precisions.append(precision)
        recalls.append(recall)

    # deal with epitopes for which we had no prediction whatsoever
    uncalled_epis = set(df['epitope']).difference(labels)
    for i in uncalled_epis:
        if not dryrun:
            recalls.append(0)
            precisions.append(0)

    recall = sum(recalls)
    precision = sum(precisions)
    fscore = sum(hmean([recalls,precisions]))

    return precision, recall, fscore


def get_clustermetrics(df, include_uncalled=True, debug_uncalled=True):
    '''Compute cluster metrics
    :param df: input data
    :type df: pandas DataFrame
    :return: scores
    :rtype: dict'''

    # Find clustered TCRs
    sub = df[~df['cluster'].isnull()]

    # Compute overall purity metrics
    stats = get_purity(sub)

    # Compute predictive metrics
    ypred = sub['cluster'].map(stats['most_frequent'])
    ytrue = sub['epitope']

    try:
        accuracy = balanced_accuracy_score(ytrue,ypred,adjusted=True)
    except ZeroDivisionError:
        print('Classification failed, zero classes')
        accuracy = 0
    precision, recall, f1score, support = precision_recall_fscore_support(ytrue, ypred, average='weighted', zero_division=0) # Total scores
    ami = adjusted_mutual_info_score(ytrue,ypred)
    rep = classification_report(ytrue, ypred, output_dict=True, zero_division=0) # Scores per epitope

    # FN are ignored above (all input TCRs are callable)!
    if include_uncalled:
        dryrun = False 
        if debug_uncalled:
            print(f"ORIGINAL recall={recall} precision={precision} f1score={f1score}")
        precision, recall, f1score = precision_recall_fscore_include_uncalled(
                df, ytrue, ypred, dryrun=dryrun, debug=debug_uncalled)

        # invalidate other values based on ytrue and ypred FIXME 
        accuracy = np.nan
        support = np.nan
        ami = np.nan
        for k in rep:
            if isinstance(rep[k], dict):
                for m in rep[k]:
                    rep[k][m] = np.nan
            else:
                rep[k] = np.nan
        if debug_uncalled:
            print(f"FIXED recall={recall} precision={precision} f1score={f1score}")

    # Compute epitope-specific metrics
    counts = {k:v for k,v in sub['epitope'].value_counts().reset_index().values.tolist()}

    epmetrics = {key:{ep:rep[ep][key] for ep in counts.keys()} for key in ['precision','recall','f1-score','support']} # Pull out epitope-specific scores from classification report
    maincluster = {ep: sub[sub['epitope']==ep]['cluster'].value_counts().index[0] for ep in sub['epitope'].unique()} # Find the largest cluster per epitope
    mosfreq=stats['most_frequent']    # Find the most frequent epitope per cluster
    clusts = {ep: '-1' if ep not in list(mosfreq.values()) else [c for c in mosfreq.keys() if mosfreq[c]==ep] for ep in counts.keys()} # Map epitopes to the clusters in which they are most frequent
    puritymap = {ep: 0 if clusts[ep]=='-1' else np.mean([len(sub[(sub['cluster']==c)&(sub['epitope']==ep)])/len(sub[sub['cluster']==c]) for c in clusts[ep]]) for ep in clusts.keys()} # Get purity per epitope
    retmap = {ep: len(df[(df['epitope']==ep)&(~df['cluster'].isnull())])/len(df[df['epitope']==ep]) for ep in counts.keys()} # Get retention scores per epitope
    consistencymap = {ep: len(sub[(sub['epitope']==ep)&(sub['cluster']==maincluster[ep])])/counts[ep] for ep in counts.keys()}  # Get consistency scores per epitope
    epmetrics['consistency']=consistencymap
    epmetrics['retention']=retmap
    epmetrics['purity']=puritymap
    
    return {'purity': np.mean(list(stats['purity_frequent'].values())),           # Purity of all clusters weighted equally (frequency)
            'purity_enriched': np.mean(list(stats['purity_enriched'].values())),  # Purity of clusters (enrichment)
            'retention': len(df[~df['cluster'].isnull()])/len(df), # Proportion of clustered TCRs
            'consistency': np.mean([(consistencymap[ep]*counts[ep])/len(sub) for ep in consistencymap.keys()]), # Proportion of an epitope assigned to a given cluster
            'ami':ami,  # Adjusted mutual information
            'accuracy':accuracy,    # Balanced accuracy
            'precision':precision,  # Precision over all epitopes
            'recall':recall,    # Recall over all epitopes
            'f1-score':f1score, # F1 over all epitopes
            'support':support,  # Support
            'mean_clustsize': np.mean(list(stats['N'].values())), # Average cluster size
            'small_clusters': stats['Nd'],    # Number of clusters with ≤ 10 members
            }, epmetrics, stats


def score(df, header):
    '''Get scores for a given set of parameters
    :param df: cluster outputs
    :type df: Pandas DataFrame
    :param header: parameters
    :type header: dict
    :return: all scores, epitope-specific scores and cluster statistics
    :rtype: Pandas DataFrame'''

    clusterscores, epscores, clusterstats = get_clustermetrics(df)  # Compute scores
    
    # Prepare dataframes
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
