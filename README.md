# tcalifornicus_projects (Oregon State University)
## This repository includes odds-and-ends associated with Tigriopus projects, in the form of scripts and instructions
### Some of these scripts have been featured in various publications during my NSF postdoctoral fellowship, including: in Molecular Ecology (https://doi.org/10.1111/mec.14973), PNAS (https://doi.org/10.1073/pnas.1819874116), and MBE (https://doi.org/10.1093/molbev/msaa008)
#### [1] TCALIF_promoter_instructions: instructions for how to extract upstream promoter regions for a genome not currently available through Bioconductor packages (based on "Bioinformatics Data Skills" by Vince Buffalo) *This is from an older version of R, now instead see bash script below*
#### [2] extract_promoters.sh: requires samtools, bamtools and gfftobed (below, edited from original version on github by RimGubaev)
#### [3] EdgeR_mRNA_data: instructions for how to perform a basic differential expression analysis in program EdgeR, with specific commands and example output for those commands
#### [4] TCALIF_mirdeep2_protocol: instructions on how to utilize the mirdeep2 pipeline for identifying miRNAs from sequencing data, with specific commands and explanations (involves mirbase mature miRNAs - mirbase_mature_jan_2018.fa)
#### [5] TCALIF_RNAseq_protocol: commands and flags for aligning mRNA reads several different ways including Salmon, Salmon + bowtie, and STAMPY + HTseq
#### [6] fasta_to_gff.pl: converts fasta file to gff3 format, attributed to Scott Cain (also available here https://github.com/scottcain/chado_test/blob/master/chado/bin/gmod_fasta2gff3.pl)
#### [7] TCALIF_TopGO: instructions for how to use TopGO in R for idenifying Gene Ontology (GO) terms associated with significantly differentially expressed genes
#### [8] edgeR_contrasts_between_times.R: instructions for how to compare different time points against each other, in combination with their controls
#### [9] mirLab: instructions for how to compare paired miRNA/mRNA-seq datasets for correlations on expression patterns
#### [10] for_loop_HMMR pipeline: instructions for how to comb through protein models for proteins containing a specific domain
#### [11] matrixEQTL.R: instructions for how to use matrixEQTL to get cis- and trans-eQTLs from snp + gene expression data
#### [12] eQTL_notes: instructions for how to format input files (snp_location and genotype) for matrixEQTL from a vcf
