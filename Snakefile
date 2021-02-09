#shell.executable("/bin/bash")
#shell.prefix("source $HOME/.bashrc; ")

import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

configfile: "config.yaml"

FULL_PATH = config['full_path']
AADIR=config['busco_dir']
OUTGROUP = config['outgroup']
IDS, = glob_wildcards(os.path.join(AADIR,"{id}.concatenated.faa"))

rule all:
	 input: 
            # MAFFT Alignment
            expand(os.path.join(AADIR,'mafft/{id}.aln',), id=IDS),
            expand(os.path.join(AADIR,'mafft/{id}.aln.trimmed',), id=IDS),
            expand(os.path.join(AADIR, 'mafft/{id}.aln.trimmed.ordered',), id = IDS),
            os.path.join(AADIR, 'mafft/concatenated_sequences.aln'),os.path.join(AADIR,'mafft/concatenated_sequences.aln.trimmed',), 
            # FASTTREE
            os.path.join(AADIR, 'trees/fasttree.nwk'),
            # RAXML (Long run time) 
            #os.path.join(FULL_PATH, AADIR, 'trees', 'raxml', 'RAxML_bestTree.EUK-MAG')    
            #os.path.join(FULL_PATH, AADIR, 'trees', 'RAxML_bestTree.EUK-MAG.rerooted.nwk') 
            expand(os.path.join(AADIR, 'mafft/{id}.aln.clustalw.out'), id = IDS)
localrules: order_add_missing,concatenated_fasta

rule mafft:
	input: aa=os.path.join(AADIR,'{id}.concatenated.faa')
	output: 
		aln=os.path.join(AADIR,'mafft/{id}.aln',)
	conda: "envs/mafft.yaml"
	shell:'''
          mafft --thread -8 --auto {input.aa} > {output.aln} 
          '''

rule trimal:
    input: os.path.join(AADIR,'mafft/{id}.aln',)
    output: os.path.join(AADIR,'mafft/{id}.aln.trimmed',)
    conda: "envs/mafft.yaml"
    shell:'''
          trimal -in {input} -out {output} -automated1
          '''

rule clustalw:
    input: os.path.join(AADIR, 'mafft/{id}.aln.trimmed')
    output: os.path.join(AADIR, 'mafft/{id}.aln.clustalw.out')
    conda: "envs/clustalo.yaml"
    shell: '''
           clustalo --infile {input} --distmat-out {output} --percent-id --full
           '''

rule order_add_missing:
    input: os.path.join(AADIR,'mafft/{id}.aln.trimmed',), os.path.join(AADIR, 'ALLMAGS.list')
    output: os.path.join(AADIR, 'mafft/{id}.aln.trimmed.ordered',)
    run: 
        all_mags = []
        with open(input[1], 'r') as f:
            for line in f:
                all_mags.append(line.strip())
        seqdict = SeqIO.to_dict(SeqIO.parse(input[0], 'fasta')) 
        # Ensure all trimmed seq lengths are the same
        lens=[]
        for i in seqdict:
            ll = len(seqdict[i].seq)
            lens.append(ll)
        assert len(set(lens))==1, 'Multiple lengths in trimmed file' 
        l=set(all_mags).difference(seqdict.keys())
        # Add empty alignments for missing sequences. 
        for i in l: 
            seqdict[i]=SeqRecord(Seq("-"*lens[0]),id=i,name=i,description=i)
        with open(output[0], "w") as outfasta:
            for record in all_mags:
                SeqIO.write(seqdict[record], outfasta,'fasta')

rule concatenated_fasta:
    input: faa= expand(os.path.join(AADIR, 'mafft/{id}.aln.trimmed.ordered',), id=IDS), mags=os.path.join(AADIR, 'ALLMAGS.list')
    output: os.path.join(AADIR, 'mafft/concatenated_sequences.aln')
    run:
        merged_seqdict = {}
        all_mags = [] 
        with open(input['mags'], 'r') as f:
            for line in f:
                all_mags.append(line.strip())
        for f in input['faa']:
            seqdict = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
            for m in all_mags:
                if m in merged_seqdict.keys():
                    append_seq = merged_seqdict[m].seq+(seqdict[m].seq)
                    merged_seqdict[m]=SeqRecord(append_seq,id=m,name=m,description=m)
                else:
                    merged_seqdict[m]=SeqRecord(seqdict[m].seq,id=m,name=m,description=m)
        with open(output[0], "w") as outfasta:
            for record in all_mags:
                SeqIO.write(merged_seqdict[record], outfasta,'fasta')

rule trimal_concat:
    input: os.path.join(AADIR,'mafft/concatenated_sequences.aln',)
    output: os.path.join(AADIR,'mafft/concatenated_sequences.aln.trimmed',)
    conda: "envs/mafft.yaml"
    shell:'''
          trimal -in {input} -out {output} -automated1
          '''

rule fast_tree: 
    input: os.path.join(AADIR,'mafft/concatenated_sequences.aln.trimmed',)
    output: os.path.join(AADIR, 'trees/fasttree.nwk')
    conda: "envs/raxmll.yaml"
    shell:'''
          FastTree -boot 100 {input} > {output}
          '''

rule raxml:
    input: os.path.join(AADIR,'mafft/concatenated_sequences.aln.trimmed',)
    output: os.path.join(FULL_PATH, AADIR, 'trees', 'raxml', 'RAxML_bestTree.EUK-MAG') 
    params: name = 'EUK-MAG', outdir = os.path.join(FULL_PATH, AADIR, 'trees', 'raxml')
    conda: "envs/raxmll.yaml"
    shell:'''
          mkdir -p {params.outdir}
          raxmlHPC-PTHREADS-SSE3 -T 16 -f a -m PROTGAMMAJTT -N 100 -n {params.name} -w {params.outdir} -s {input} -p 42 -x 42
          '''

rule re_root_raxml:
    input: os.path.join(FULL_PATH, AADIR, 'trees', 'raxml', 'RAxML_bestTree.EUK-MAG')
    output: os.path.join(FULL_PATH, AADIR, 'trees', 'RAxML_bestTree.EUK-MAG.rerooted.nwk')
    conda: "envs/raxmll.yaml"
    params: outgroup = OUTGROUP
    shell:'''
          nw_reroot {input} {params.outgroup} > {output}
          '''
