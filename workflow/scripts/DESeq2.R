# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input


#############################################################
#----------------- Loading the librairies -------------------
#############################################################
library(readr)
library(DESeq2)

# Variables coming from the snakemake
output_dir <- snakemake@params[["out_dir"]]
design_file <- snakemake@input[["samples"]]
count_file <- snakemake@input[["quant"]]
comparison_file <- snakemake@input[["comparisons"]]

# create the directory that will contain the results
dir.create(output_dir, showWarnings=FALSE)


#############################################################
#------------- Importing data and information ---------------
#############################################################
# Loading samples information from the design file.
sampleTable <- read.table(
    design_file, header=TRUE,
    row.names="sample", check.names=FALSE
)
conditions <- unique(sampleTable$condition)
samples <- rownames(sampleTable)
sampleTable$condition <- factor(sampleTable$condition, levels=conditions)

# Loading the comparisons
comparisons <- read.table(
    comparison_file, header=TRUE
)

cts <- read.csv(count_file,sep="\t",row.names="gene_id")


#############################################################
#--------------- Creating the DESeq2 object -----------------
#############################################################
# Create the DESeq2 object with all the conditions
dds <- DESeqDataSetFromMatrix(
    cts,
    sampleTable,
    ~condition
)

# pre filtering the count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# put the levels for each columns
dds$condition <- factor(dds$condition, conditions)

# call the function for differential gene expression analysis
dds <- DESeq(dds)


#############################################################
#-------------- Looping over all comparisons ----------------
#############################################################
for (row in 1:nrow(comparisons)) {

    condition1 <- comparisons[row, "cdn1"]
    condition2  <- comparisons[row, "cdn2"]
    exp <- sprintf("%s-%s", condition1, condition2)

    res <- results(
        dds,
        contrast = c(
            "condition",
            condition1,
            condition2
        ),
        #pAdjustMethod = "fdr"
    )

    # transform the result to data frame
    res_df <- as.data.frame(res)

    # Writing results to file
    fname <- paste(
        output_dir,
        paste(exp, "csv", sep='.'),
        sep='/'
    )

    write.csv(
        res_df,
        file=fname,
        quote=FALSE,
    )
}