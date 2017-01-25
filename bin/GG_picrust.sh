#!/bin/bash

echo "Green Genes Data and Picrust"

pick_closed_reference_otus.py -i IMG.fasta -o IMG -p otu_picking_params_97.txt -r gg_13_5_otus/rep_set/97_otus.fasta -t gg_13_5_otus/taxonomy/97_otu_taxonomy.txt

predict_metagenomes.py -i IMG/otu_table.biom -o IMG/metagenome_predicted.biom

biom convert -i metagenome_predicted.biom -o metagenome_predicted.txt --to-tsv

biom convert -i otu_table.biom -o otu_table.txt --to-tsv

categorize_by_function.py -i metagenome_predicted.biom -c KEGG_Pathways -l 3 -o predicted_patyways.L3.biom

biom convert -i predicted_patyways.L3.biom -o predicted_patyways.L3.txt --to-tsv
