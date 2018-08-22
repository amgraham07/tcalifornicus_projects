library(edgeR)

#import the expression data
x <- read.delim("hypo.txt",row.names="Gene")

#import the data frame
targets <- read.delim("targets.txt",row.names="Sample")

#create a group by concatenating the two columns in the data frame:
Group <- factor(paste(targets$Treat,targets$Time,sep="."))
cbind(targets,Group=Group)

#create the usual edgeR object:
y <- DGEList(counts=x,group=Group)
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+Group) #make sure the ~0+ is found in this case of multiple groups.
colnames(design) <- levels(Group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

#make a list of contrasts:
my.contrasts <- makeContrasts( Anox24hvs3h = (anox.24h-normo.24h)-(anox.3h-normo.3h),levels=design)

#run the model
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"Anox24hvs3h"])

#check the top genes
topTags(qlf)

#calculate FDR values for the Pvalues, then combine them to the table:
fdr=p.adjust(qlf$table$PValue, method="fdr")
new.output=cbind(qlf$table, fdr)

#export the final table
write.table(new.output, "Anox24hvs3h_output.txt")
