#!/bin/bash
set -euo pipefail  

# Input files
ARAPORT_RAW="Araport11_pep_20220914"
ARAPORT_CLEAN="Araport11_pep_20220914.fa"
BIDENTIS_QUERY="F_bidentis.proteins.fa"
GYNANDROPSIS_QUERY="gyn-vHIC-CDS-prot-new.fasta.fa"
GFF3_FILE="gene.mRNA.UTR.Fbidentis_MAKER_Helixer.AGAT.LongIso.gff3"

# Step 1: Clean Araport peptide file (remove non-ASCII characters)
echo "Cleaning Araport peptide file..."
perl -pe 's/[^\x00-\x7F]+/ /g' "$ARAPORT_RAW" > "$ARAPORT_CLEAN"

# Step 2: Make BLAST database
echo "Creating BLAST database..."
makeblastdb -in "$ARAPORT_CLEAN" -dbtype prot -out "$ARAPORT_CLEAN" -parse_seqids

# Step 3: Run BLASTP (Bidentis vs Arabidopsis)
echo "Running BLASTP for Bidentis..."
blastp -db "$ARAPORT_CLEAN" -query "$BIDENTIS_QUERY" \
    -out bidentis_arabidopsis.out -outfmt 7 -evalue 1e-5 \
    -max_target_seqs 5 -word_size 3 -num_threads 3

# Step 4: Run BLASTP (Gynandropsis vs Arabidopsis)
echo "Running BLASTP for Gynandropsis..."
blastp -db "$ARAPORT_CLEAN" -query "$GYNANDROPSIS_QUERY" \
    -out gyn_to_arabidopsis.out -outfmt 7 -evalue 1e-5 \
    -max_target_seqs 5 -word_size 3 -num_threads 3

# Step 5: Get best hits
echo "Extracting best hits..."
python get_besthit.py bidentis_arabidopsis.out besthit_bidentis_arabidopsis.out
python get_besthit.py gyn_to_arabidopsis.out besthit_gyn_to_arabidopsis.out

# Step 6: Change names for Bidentis
echo "Changing gene names for Bidentis best hits..."
python bid_change_name.py "$GFF3_FILE" besthit_bidentis_arabidopsis.out new_names_besthit_bidentis_arabidopsis.out

echo "All steps completed successfully!"
