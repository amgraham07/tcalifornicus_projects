# tcalifornicus_projects (OSU)
## This repository includes odds-and-ends associated with Tigriopus projects, in the form of scripts and instructions
### Some of these scripts have been featured in various publications during my NSF postdoctoral fellowship, including: in Molecular Ecology (https://doi.org/10.1111/mec.14973), PNAS (https://doi.org/10.1073/pnas.1819874116)
#### [1] TCALIF_promoter_instructions: instructions for how to extract upstream promoter regions for a genome not currently available through Bioconductor packages (based on "Bioinformatics Data Skills" by Vince Buffalo) *This is from an older version of R, now instead see bash script below*
#### [2] extract_promoters.sh: requires samtools and bamtools, edited from original version on github by RimGubaev
#### [3] EdgeR_mRNA_data: instructions for how to perform a basic differential expression analysis in program EdgeR, with specific commands and example output for those commands
#### [4] TCALIF_mirdeep2_protocol: instructions on how to utilize the mirdeep2 pipeline for identifying miRNAs from sequencing data, with specific commands and explanations (involves mirbase mature miRNAs - mirbase_mature_jan_2018.fa)
#### [5] TCALIF_RNAseq_protocol: commands and flags for aligning mRNA reads to a reference using Bowtie2, or quasi-align the reads using SALMON, plus quantification using SALMON
#### [6] TCALIF_TopGO: instructions for how to use TopGO in R for idenifying Gene Ontology (GO) terms associated with significantly differentially expressed genes
#### [7] edgeR_contrasts_between_times.R: instructions for how to compare different time points against each other, in combination with their controls
#### [8] mirLab: instructions for how to compare paired miRNA/mRNA-seq datasets for correlations on expression patterns
#### [9] for_loop_HMMR pipeline: instructions for how to comb through protein models for proteins containing a specific domain
