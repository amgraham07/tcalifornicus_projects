# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');

## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name <- read.delim("genotypes.txt",row.names="Name");
snps_location_file_name <- read.delim("snp_locations.txt",row.names="snp");

# Gene expression file name
#expression_file_name <- read.delim("eQTL_logcpm_control.txt",row.names="Name");
expression_file_name <- read.delim("eQTL_logcpm_heat.txt",row.names="Name");
gene_location_file_name <- read.delim("gene_locations.txt",row.names = "geneid");

# Covariates file name
#Set to character() for no covariates
covariates_file_name <- read.delim("covariate.txt",row.names="Name");
#covariates_file_name <- character();


# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-3;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e5;

## Load genotype data ## slicing is for memory purposes
snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in slices of 2,000 rows
snps$LoadFile("genotypes.txt");

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t"; # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in slices of 2,000 rows
#gene$LoadFile("eQTL_logcpm_control.txt");
gene$LoadFile("eQTL_logcpm_heat.txt");

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile("covariate.txt");
}

## Run the analysis
snpspos = read.table("snp_locations.txt", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("gene_locations.txt", header = TRUE, stringsAsFactors = FALSE);
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Make the histogram of local and distant p-values
plot(me)

## write output
write.table(me$cis$eqtls,"cis_p2_co_heat_1e5.txt")
#write.table(me$trans$eqtls,"trans_p2_co_heat.txt")

#CIS
## install plyr, identify which snps are sig associated with which genes at a max FDR of 10%
## print table of only the lead SNP per gene, also add the MAF for every SNP in table
cis_eQTL_res = me$cis$eqtls
cis_eQTL_res = cis_eQTL_res[cis_eQTL_res$FDR < 0.1,]

top_eqtls = cis_eQTL_res[order(cis_eQTL_res$pvalue),]
top_eqtls = top_eqtls[!duplicated(top_eqtls$gene),]

mafs = apply(as.matrix(SNP_file_name),1,mean)/2
mafs = data.frame(snps=names(mafs), maf = mafs)

top_eqtls = merge(top_eqtls, mafs, by="snps")
top_eqtls = top_eqtls[order(top_eqtls$FDR),]

head(top_eqtls)
write.table(top_eqtls,"top_cis_p2_co_heat_1e5.txt")

#TRANS
## install plyr, identify which snps are sig associated with which genes at a max FDR of 10%
## print table of only the lead SNP per gene, also add the MAF for every SNP in table
trans_eQTL_res = me$trans$eqtls
trans_eQTL_res = trans_eQTL_res[trans_eQTL_res$FDR < 0.1,]

trans_top_eqtls = trans_eQTL_res[order(trans_eQTL_res$pvalue),]
trans_top_eqtls = trans_top_eqtls[!duplicated(trans_top_eqtls$gene),]

head(trans_top_eqtls)
write.table(trans_top_eqtls,"top_trans_p2_co_heat.txt")
