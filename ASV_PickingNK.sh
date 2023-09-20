#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=245G
#SBATCH --error=/lustre/work/benson/nkorth/SAP/ASVPicking_S766/job.%J.err
#SBATCH --output=/lustre/work/benson/nkorth/SAP/ASVPicking_S766/job.%J.out
#SBATCH --job-name=step1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nate.korth@gmail.com
#SBATCH --partition=benson
#SBATCH --nodes=4

#cd /lustre/work/benson/nkorth/SAP/ASVPicking_S766/Data
#gunzip *.gz

cd /lustre/work/benson/nkorth/SAP/ASVPicking_S766
#wget -O 'silva-132-99-515-806-nb-classifier.qza' 'https://data.qiime2.org/2019.4/common/silva-132-99-515-806-nb-classifier.qza'

module load qiime2/2019.1
#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path SAP_S1_manifest.csv --output-path SAP_S1_paired-end-demux.qza --input-format PairedEndFastqManifestPhred33
#qiime tools peek SAP_S1_paired-end-demux.qza
#qiime demux summarize --i-data SAP_S1_paired-end-demux.qza --o-visualization SAP_S1_paired-end-demux.qzv
#qiime dada2 denoise-paired --i-demultiplexed-seqs SAP_S1_paired-end-demux.qza --o-table SAP_S1_table --o-representative-sequences SAP_S1_rep-seqs --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 220 --p-trunc-len-r 160 --p-hashed-feature-ids --o-denoising-stats stats-dada2.qza --p-n-threads 0 --verbose
#qiime feature-table summarize --i-table SAP_S1_table.qza --o-visualization SAP_S1_table.qzv --m-sample-metadata-file SAP_S1_metadata.tsv
#qiime feature-table tabulate-seqs --i-data SAP_S1_rep-seqs.qza --o-visualization SAP_S1_rep-seqs.qzv
#qiime alignment mafft --i-sequences SAP_S1_rep-seqs.qza --o-alignment SAP_S1_aligned-rep-seqs.qza --verbose --p-n-threads -0
#qiime alignment mask --i-alignment SAP_S1_aligned-rep-seqs.qza --o-masked-alignment SAP_S1_masked-aligned-rep-seqs.qza --verbose
#qiime phylogeny fasttree --i-alignment SAP_S1_masked-aligned-rep-seqs.qza --o-tree SAP_S1_unrooted-tree.qza --verbose --p-n-threads 1
#qiime phylogeny midpoint-root --i-tree SAP_S1_unrooted-tree.qza --o-rooted-tree SAP_S1_rooted-tree.qza
#qiime feature-classifier classify-sklearn --i-classifier silva-132-99-515-806-nb-classifier.qza --i-reads SAP_S1_rep-seqs.qza --o-classification SAP_S1_silva-taxonomy.qza --p-n-jobs -2 --verbose
#qiime metadata tabulate --m-input-file SAP_S1_silva-taxonomy.qza --o-visualization SAP_S1_silva-taxonomy.qzv
#qiime taxa barplot --i-table SAP_S1_table.qza --i-taxonomy SAP_S1_silva-taxonomy.qza --m-metadata-file SAP_S1_metadata.tsv --o-visualization SAP_S1_silva-taxa-bar-plots.qzv
#qiime tools export --input-path SAP_S1_table.qza --output-path exported
#mv ./exported/feature-table.biom ./exported/SAP_S1_OTU_Table.biom
#qiime tools export --input-path SAP_S1_silva-taxonomy.qza --output-path exported
#mv ./exported/taxonomy.tsv ./exported/SAP_S1_silva_taxonomy.tsv
#qiime tools export --input-path SAP_S1_rep-seqs.qza --output-path exported
#mv ./exported/dna-sequences.fasta ./exported/SAP_S1_rep_seqs.fasta
#qiime tools export --input-path SAP_S1_rooted-tree.qza --output-path exported
#mv ./exported/tree.nwk ./exported/SAP_S1_tree.nwk
#cd ./exported
#cp SAP_S1_silva_taxonomy.tsv SAP_S1_silva_taxonomy_biom.tsv
#sed -i '1d' SAP_S1_silva_taxonomy_biom.tsv
#sed -i '1i #OTUID taxonomy confidence' SAP_S1_silva_taxonomy_biom.tsv
#biom add-metadata -i SAP_S1_OTU_Table.biom -o SAP_S1_OTU_Table_silva_taxonomy.biom --observation-metadata-fp SAP_S1_silva_taxonomy_biom.tsv --sc-separated taxonomy
#biom convert -i SAP_S1_OTU_Table_silva_taxonomy.biom -o SAP_S1_OTU_Table_silva_taxonomy.tsv --to-tsv --header-key taxonomy --output-metadata-id 'Consensus Lineage'

#cd ../
#Samples must have more than 1500 reads
#qiime feature-table filter-samples  --i-table SAP_S1_table.qza  --p-min-frequency 1500 --o-filtered-table temp.qza
#ASV must be present in more than 2 samples
#qiime feature-table filter-features  --i-table temp.qza  --p-min-samples 3  --o-filtered-table temp2.qza
#ASV must have more than 15 total reads
#qiime feature-table filter-features  --i-table temp2.qza  --p-min-frequency 15  --o-filtered-table temp3.qza

#qiime taxa collapse --i-table temp3.qza --i-taxonomy SAP_S1_silva-taxonomy.qza --p-level 6 --o-collapsed-table SAP_S1_genus_table.qza

#qiime tools export --input-path SAP_S1_genus_table.qza --output-path exported
#mv ./exported/feature-table.biom ./exported/SAP_S1_OTU_genus_table.biom

cd exported
#biom add-metadata -i SAP_S1_OTU_genus_table.biom -o SAP_S1_OTU_genus_Table_silva_taxonomy.biom --observation-metadata-fp SAP_S1_silva_taxonomy_biom.tsv --sc-separated taxonomy
#biom convert -i SAP_S1_OTU_genus_Table_silva_taxonomy.biom -o SAP_S1_OTU_genus_Table_silva_taxonomy.tsv --to-tsv --header-key taxonomy --output-metadata-id 'Consensus Lineage'
biom convert -i SAP_S1_OTU_genus_Table_silva_taxonomy.tsv -o SAP_S1_OTU_genus_Table_silva_taxonomy.json.biom --table-type="OTU table" --to-json

#cd ..

#qiime tools export --input-path temp3.qza --output-path exported
#mv ./exported/feature-table.biom ./exported/SAP_S1_OTU_trim_Table.biom

#cd exported

#biom add-metadata -i SAP_S1_OTU_trim_Table.biom -o SAP_S1_OTU_trim_Table_silva_taxonomy.biom --observation-metadata-fp SAP_S1_silva_taxonomy_biom.tsv --sc-separated taxonomy
#biom convert -i SAP_S1_OTU_trim_Table_silva_taxonomy.biom -o SAP_S1_OTU_trim_Table_silva_taxonomy.tsv --to-tsv --header-key taxonomy --output-metadata-id 'Consensus Lineage'


#biom convert -i SAP_S1_OTU_trim_Table_silva_taxonomy.tsv -o SAP_S1_OTU_trim_Table_silva_taxonomy.json.biom --table-type="OTU table" --to-json













