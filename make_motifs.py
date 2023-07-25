import os, sys,subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
from collections import Counter

def get_motif(sequences,name,output_folder):
    instances = [Seq(x) for x in sequences]
    record = [SeqRecord(seq, id=str(s), description='cluster_%s_instance_%s'%(name,str(s)),
                        annotations={"molecule_type": "protein"}) for s, seq in enumerate(instances)]

    if not 'Motifs' in os.listdir(output_folder):
        os.mkdir(os.path.join(output_folder,'Motifs'))
    
    dest = '%s/Motifs/%s.fa'%(output_folder,name)
    SeqIO.write(record,dest,'fasta')
    out_file= '%s/Motifs/%s_msa.fa'%(output_folder,name)

    # cmd='muscle -super5 {} -output {}'.format(dest,out_file) # MUSCLE 5.1.osxarm64
    cmd='muscle -in {} -out {}'.format(dest,out_file) # MUSCLE v3.8.31

    if sys.platform.lower() == 'darwin':
        subprocess.call(cmd, shell=True, executable='/bin/zsh')
    else:
        subprocess.call(cmd, shell=True, executable='/bin/bash')

    out_logo = '%s/%s.eps'%(output_folder,name)
    cmd = 'weblogo -f {} -o {} -F eps -P {} -s large'.format(out_file, out_logo,name)
    os.system(cmd)
    os.system('rm {}'.format('%s/Motifs/*.fa'%(output_folder)))


if __name__=="__main__":
    
    models=['clustcr',
            'GIANA',
            'gliph2',
            'hamming',
            'ismart',
            'tcrdist3',
            'length',
            'random'
            ]
    names= ['ClusTCR',
            'GIANA',
            'GLIPH2',
            'Hamming',
            'iSMART',
            'tcrdist3',
            'Length',
            'Random']
    n=250
    root= 'results/Motifs/20230720_081643'
    for m, model in enumerate(models):

        print(model)
        data = pd.read_csv(os.path.join(root,'%s_beta.csv'%(model)))
        topn = data['cluster'].value_counts().index[:n]

        # lengths =[[len(cdr3) for cdr3 in data[data['cluster']==c]['cdr3.beta'].values] for c in topn]
        for i, c in enumerate(topn):
            sub=data[data['cluster']==c]
            sub['length']= [len(cdr3) for cdr3 in sub['cdr3.beta'].values]
            mode=Counter(sub['length']).most_common(1)[0][0]
            ep=Counter(sub['epitope']).most_common(1)[0][0]
            if ep not in os.listdir('results/Motifs/20230720_081643/'):
                os.mkdir('results/Motifs/20230720_081643/%s'%ep)
            sub2= sub[sub['length']==mode]
            seqs = sub2['cdr3.beta'].dropna().values
            get_motif(seqs,'%s_C%s_%s'%(model,str(i+1),ep),os.path.join(os.getcwd(),'results/Motifs/20230720_081643/%s'%(ep)))

