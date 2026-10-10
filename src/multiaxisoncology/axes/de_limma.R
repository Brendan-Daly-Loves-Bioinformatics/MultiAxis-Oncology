# Differential expression: tumor vs normal (Toil-harmonized TPM), config-driven
# Usage: Rscript de_limma.R <expr_csv.gz> <probemap> <out_csv>
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
expr_path <- args[1]; pm_path <- args[2]; out_csv <- args[3]

expr <- read.csv(gzfile(expr_path), stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
mat <- as.matrix(expr)

pm <- read.delim(pm_path, stringsAsFactors = FALSE)
map <- pm$gene[match(rownames(mat), pm$id)]
keep <- !is.na(map)
mat <- mat[keep, ]; map <- map[keep]

# Collapse duplicate gene symbols by max mean expression
mn <- rowMeans(mat)
ord <- order(map, -mn); mat <- mat[ord, ]; map <- map[ord]
d <- !duplicated(map); mat <- mat[d, ]; rownames(mat) <- map[d]
cat("Genes after symbol collapse:", nrow(mat), "\n")

samples <- colnames(mat)
group <- ifelse(grepl("-11$", samples) | grepl("^GTEX", samples), "normal",
         ifelse(grepl("-0[16]$", samples), "tumor", NA))
table(group, useNA = "ifany")
keep_s <- !is.na(group)
mat <- mat[, keep_s, drop = FALSE]; group <- group[keep_s]

expr_ok <- rowSums(mat > 1) >= 0.1 * ncol(mat)
cat("Genes expressed filter:", sum(expr_ok), "\n")
mat <- mat[expr_ok, ]

design <- model.matrix(~ 0 + group)
colnames(design) <- c("normal", "tumor")
fit <- lmFit(mat, design)
cont <- makeContrasts(tumor - normal, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))

tt <- topTable(fit2, number = Inf, adjust.method = "BH")
tt$gene <- rownames(tt)
tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(tt, out_csv, row.names = FALSE)
cat("Up in tumor (FDR<0.05, logFC>1):", sum(tt$adj.P.Val < 0.05 & tt$logFC > 1), "\n")
cat("Down in tumor (FDR<0.05, logFC<-1):", sum(tt$adj.P.Val < 0.05 & tt$logFC < -1), "\n")
