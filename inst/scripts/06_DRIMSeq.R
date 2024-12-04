#' ---
#' title: "Differential transcript usage using DRIMSeq"
#' author: "Nicolas Delhomme for BN Bioinformatics"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Rationale
#' ## Selecting tools
#' DTU is still a very young topic in expression profiling, despite some tools having been available for a decade.
#' 
#' It is hard to know what tool to use, even if there are not that many available yet.
#' 
#' Figuring it out is a good example of "technological watch" in bioinformatics. Selecting a (few) tool(s) will involve:
#' 
#' 1. Assessing its (their) relevance through a literature search, possibly identifying some benchmark report.
#' 2. Asserting the code availability, and its quality (is it maintained and actively developed, or was it a one of for publication)
#' 3. Installing the tool (that's an easy filter, if a tool does not install easily, that's a no-no)
#' 4. Looking up the tool documentation (this could be part of 2. and could be swapped with 3., but 3. is commonly fast). If the documentation is scarce or hard to find, that's minus points.
#' 
#' Eventually, you will probably end up with a few tools that fulfill all criteria above and a lot of the selection process has to be
#' pragmatic. Even if a tool is not the "best" according to benchmarks, but is easy to install, run and has good doc, it's probably the better choice (for the sake of your time).
#' 
#' ## A practical example
#' 
#' There is one good benchmark published by [Love _et al._](https://f1000research.com/articles/7-952/v3). 
#' 
#' For DTU, they compare four tools:
#' 1. `DEXSeq` 
#' 2. `DRIMSeq`
#' 3. `RATs`
#' 4. `SUPPA2`
#' 
#' `IsoformSwitchAnalyzeR` is mentioned in the review, but not assessed.
#' 
#' For DGE and DTE, they compared:
#' 
#' 1. `DESeq2` 
#' 2. `EBSeq`
#' 3. `edgeR` (multiple flavours)
#' 4. `limma`
#' 5. `SAMseq`
#' 6. `sleuth`
#' 
#' ### Tool selection
#' 
#' For the DTU, the benchmark conclusions are:
#' 
#' **RATs and SUPPA2 control their FDR well, but have a lower True Positive Rate than DEXSeq and DRIMSeq.**
#' 
#' Filtering the input data (signal-to-noise ratio) does help bringing the FDR closer to its target.
#' 
#' As an answer to point 1. above, it seems all tools can be considered at that stage. 
#' 
#' Moving to point 2. all four tools have repositories that are easy to access. All, but SUPPA2, are `R` packages. SUPPA2 is written in `python`. There is nothing here that would immediately red flag a tool.
#' 
#' With regards to point 3. neither, they are all easy to install and run.
#' 
#' As for point 4. DEXSeq, and DRIMSeq have good documentation. RATs, too, but it somewhat requires more time to patch it together and it is fairly succinct. There are some instructions to import data but it is more or less restricted to importing `kallisto`.
#' 
#' Where does that leave us? **Note** that this is the most **pragmatic** part of the process, it does not imply anything about the tools' performance.
#'
#' 1. `DEXSeq` - not a DTU tool per-se
#' 2. `DRIMSeq` - a candidate that fits well our pipeline. Avail in R and documented.
#' 3. `RATs` - the documentation is limited and would require effort to figure out how to import our `salmon` data.
#' 4. `SUPPA2` - the tool is written in `python` so less fitted to our environment.
#' 
#' There is also `IsoformSwitchAnalyzeR`, as a well documented R package, but we lack an in-depth comparison to other tools.
#' 
#' As a final note, the benchmark manuscript has an associated [GitHub repository](https://github.com/mikelove/swimdown) which contains all the code to actually run `SUPPA2` and import the results in `R`. Similarly, it contioans the code to use `RATs` on `salmon` output. Plus all the code necessary to reproduce the whole benchmarking analysis, so definitely a hugely useful resource!
#' 
#'  ## Decision
#' 
#'  Let's give a shot at `DRIMSeq`. To retrieve the help page, we only need to run
#'  `vignette(package="DRIMSeq")` to find the available documentation in that package, before accessing it: `vignette("DRIMSeq",package="DRIMSeq")`
#' 
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(DRIMSeq)
  library(PasillaTranscriptExpr)
})

#' # Data
data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")

#' * Load metadata
(pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"),
                               header = TRUE, as.is = TRUE))

#' * Load counts
pasilla_counts <- read.table(file.path(data_dir, "counts.txt"),
                             header = TRUE, as.is = TRUE)

dim(pasilla_counts)
str(pasilla_counts)
dplyr::glimpse(pasilla_counts)

#' * `dmDSdata`
pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName,
                              group = pasilla_metadata$condition)

(d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples))

head(counts(d), 3)
head(DRIMSeq::samples(d), 3)

#' * A visual look at the number of transcripts per gene
plotData(d)

#' How many genes have more than one transcript?
(tab <- table(table(counts(d)$gene_id)))
sum(tab[-1])

#' * Subset
#' Here we just take a subset for the sake of the tutorial
gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
length(gene_id_subset)
(d <- d[names(d) %in% gene_id_subset, ])

#' # Filtering
#' As reported by [Love et al.](https://f1000research.com/articles/7-952/v3) and implemented in `DRIMSeq`, filtering for low expression (_i.e._ signal-to-noise) is an essential step.
#' 
#' There are 4 filters:
#' 
#' 1. `min_samps_gene_expr`
#' 2. `min_samps_feature_expr`
#' 3. `min_gene_expr`
#' 4. `min_feature_expr`
#' 
(d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
               min_gene_expr = 10, min_feature_expr = 10))

#' # Precision estimation
#' To calculate the precision estimate, we face the same problem as for the Differential Gene Expression (DGE), _i.e._ that we have few samples.
#'
#' However as for the DGE, we have many genes so we could use that information to calculate a common dispersion, assuming the precision is a constant for almost all genes.
#' 
#' Unlike for gene expression, this is likely to be too strong of an assumption for DTU. Rather the authors of DRIMSeq, decided on a moderated approach.
#' 
#' In that approach, 10% of the genes (by default) are selected at random to estimate the common dispersion.
#' 
#' Prior to estimating the precision, we need to define the design matrix. 
(design_full <- model.matrix(~ group, data = DRIMSeq::samples(d)))

#' Then we can calculate the precision. This is the time consuming step of DRIMSeq, which can be parallelised 
#' on a real dataset (`BPPARAM = BiocParallel::MulticoreParam()`).
(d <- dmPrecision(d, design = design_full))

#' We can take a look at the content
head(mean_expression(d), 3)
common_precision(d)
head(genewise_precision(d),3)

plotPrecision(d)

#' # Proportion estimation
#' _i.e._ fitting the model
(d <- dmFit(d, design = design_full, verbose = 1))

#' Taking a look at the proportion
head(proportions(d))

#' And the the Dirichlet-multinomial coefficients (at the gene level)
head(coefficients(d))

#' Also looking at the beta-binomial transcript results (_a.k.a._ "feature")
head(coefficients(d), level = "feature")

#' # Testing for DTU
#' This is done by fitting two models, a `null` and a `full`. Then likelihood ratio statistics are used to test for differences between the two models.
d <- dmTest(d, coef = "groupKD", verbose = 1)

#' ## Results
#' * gene level
head(results(d), 3)

#' * transcript level
head(results(d, level = "feature"), 6)

#' ## Visualisation
#' * Assessment
plotPValues(d)
plotPValues(d,level="feature")

#' * Results
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]
top_gene_id <- res$gene_id[1]
plotProportions(d, gene_id = top_gene_id, group_variable = "group")

plotProportions(d, gene_id = top_gene_id, group_variable = "group",
                plot_type = "lineplot")

plotProportions(d, gene_id = top_gene_id, group_variable = "group",
                plot_type = "ribbonplot")

#' # Two-stage test
#' One can use the `stageR` package to set up a two stage test. During the first step, a screening test, one will select genes of interest.
#' Then in the confirmation test, for each gene separately, the transcript p-values will be adjusted.
#' 
suppressPackageStartupMessages(library(stageR))

#' Assign gene-level pvalues to the screening stage
pScreen <- results(d)$pvalue
names(pScreen) <- results(d)$gene_id

#' Assign transcript-level p-values to the confirmation stage
pConfirmation <- matrix(results(d, level = "feature")$pvalue, ncol = 1)
rownames(pConfirmation) <- results(d, level = "feature")$feature_id

#' Create the gene-transcript mapping
tx2gene <- results(d, level = "feature")[, c("feature_id", "gene_id")]

#' Create the stageRTx object and perform the stage-wise analysis
stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.05)

#' Taking a look
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)
padj <- getAdjustedPValues(stageRObj, order = TRUE,
                           onlySignificantGenes = FALSE)
head(padj)

#' # Conclusion
#' 
#' And that's it, we have unlocked DGE, DTE and DTU! Cool biology and data science awaits us :-)
#' 
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
