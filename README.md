# Fish-microbiome
Gut microbial composition of the cyprinids Cyprinella lutrensis (Red shiner) and Notropis stramineus (Sand shiner): insights from wild fish populations

## List of Files:  
* FishMicrobiomeAnalysis_dada2_rscript.R: The R script submitted to the Univeristy of Nebraska Lincoln's (UNL) Holland Commputing Center (HCC) cluster to run the DADA2 pipeline automatically.
* mothur_output.txt: Mothur was used to create a phylogenetic tree of the ASV sequences in an interactive session on UNL's HCC clusters. This text file contains a record of the commands and output generated in that session.
* FishMicrobeAnalysis.Rmd: A markdown file focusing on combining the data of three sequencing runs into a single data set and assessing the negative controls of each PCR 96 well plate.
* FishMicrobeAnalysis.html: The rendered html file of FishMicrobeAnalysis.Rmd.
* FishMicrobeAnalysis3.qmd: A quarto document focused specifically on data that will be used in the final analysis. The ASV table, taxonomy table, and meta data are combined into a phyloseq object and quality control filtering is accomplished.
* FishMicrobeAnalysis3.html: The rendered html file of FishMicrobeAnalysis3.qmd
* ks_FishMicrobe_diversity_version3.qmd: A quarto document in which most statistical analysis are completed, including diversity metrics, environemental association tests, differential abundance tests, and more.
* ks_FishMicrobe_diversity_version3.html: The rendered html file of ks_FishMicrobe_diversity_version3.qmd
