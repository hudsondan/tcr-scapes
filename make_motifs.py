import os, sys,subprocess, argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
from collections import Counter

def get_motif(sequences, name, output_folder):
    '''Produce cluster motifs from input results folder
    :param sequences: input TCR sequences
    :type sequences: list
    :param name: cluster, model and epitope identifier
    :type name: str
    :param output_folder: directory for msa and motif files

    '''
    instances = [Seq(x) for x in sequences] # Create biopython Seq instances
    record = [SeqRecord(seq, id=str(s), description='cluster_%s_instance_%s'%(name,str(s)),
                        annotations={"molecule_type": "protein"}) for s, seq in enumerate(instances)]  # Convert to SeqRecord 
    
    # Initialise motifs folder
    if not 'Motifs' in os.listdir(output_folder):
        os.mkdir(os.path.join(output_folder,'Motifs'))
    
    dest = '%s/Motifs/%s.fa'%(output_folder,name)
    SeqIO.write(record,dest,'fasta')    # Write FASTA file for MUSCLE
    out_file= '%s/Motifs/%s_msa.fa'%(output_folder,name)


    cmd='muscle -super5 {} -output {}'.format(dest,out_file) # MUSCLE 5.1.osxarm64
    # cmd='muscle -in {} -out {}'.format(dest,out_file) # MUSCLE v3.8.31

    if sys.platform.lower() == 'darwin':
        subprocess.call(cmd, shell=True, executable='/bin/zsh')
    else:
        subprocess.call(cmd, shell=True, executable='/bin/bash')

    out_logo = '%s/%s.eps'%(output_folder,name)
    cmd = 'weblogo -f {} -o {} -F eps -P {} -s large'.format(out_file, out_logo,name)
    os.system(cmd)
    os.system('rm {}'.format('%s/Motifs/*.fa'%(output_folder)))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='ClustOx: Produce CDR3 motifs for clustered instances')
    parser.add_argument('-i', '--input', required=False, type=str, default='results_publication/Motifs/20230720_081643',
                        help='Set the input cluster folder, for example results/YYYMMDD_HHMMSS')
    parser.add_argument('-n', '--n_clusters', required=False, type=int, default=20,
                        help='Set the number of clusters for which logos will be produced')
    parser.add_argument('-cs', '--chain_selection', required=False, type=str, default='beta',
                        help='Set the chain selection for which logos will be produced')
    
    models=['clustcr',
            'GIANA',
            'gliph2',
            'hamming',
            'ismart',
            'tcrdist3',
            'vcluster',
            'length',
            'random'
            ]
    names= ['ClusTCR',
            'GIANA',
            'GLIPH2',
            'Hamming',
            'iSMART',
            'tcrdist3',
            'vcluster',
            'Length',
            'Random']
    
    args=parser.parse_args()
    root= args.input
    n= args.n_clusters
    chain_selection = args.chain_selection

    if chain_selection not in ['alpha','beta']:
        raise ValueError('Select a chain from ["alpha","beta"]')

    for m, model in enumerate(models):
        data = pd.read_csv(os.path.join(root,'%s_%s.csv'%(model,chain_selection)))
        topn = data['cluster'].value_counts().index[:n]

        for i, c in enumerate(topn):
            sub=data[data['cluster']==c]
            sub['length']= [len(cdr3) for cdr3 in sub['cdr3.%s'%(chain_selection)].values]
            mode=Counter(sub['length']).most_common(1)[0][0]
            ep=Counter(sub['epitope']).most_common(1)[0][0]
            if ep not in os.listdir(root):
                os.mkdir('%s/%s'%(root,ep))
            sub2= sub[sub['length']==mode]
            seqs = sub2['cdr3.%s'%(chain_selection)].dropna().values
            get_motif(seqs,'%s_C%s_%s_%s'%(model,str(i+1),ep,chain_selection),os.path.join(os.getcwd(),'%s/%s'%(root,ep)))

