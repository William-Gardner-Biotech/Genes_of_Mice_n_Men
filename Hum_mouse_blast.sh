#!/bin/bash

module load bio

# Database path of your human proteins
# database=~/YOURDATABASE

makeblastdb -dbtype 'prot' -in $database -out human_proteom.fasta

blastp -query mouse_proteome.fasta -db human_proteome.fasta \
    -out mouse_vs_human_blast_results.tab \
    -evalue 1e-20 \
    -outfmt 6\
    -max_target_seqs 1
