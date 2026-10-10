# DepMap 22Q2 CRISPR dependency: lineage-selective dependencies (config-driven)
# Usage: Rscript depmap_dependency.R <crispr_rds> <metadata_rds> <lineage> <subtype> <sub_subtype> <out_csv>
# subtype / sub_subtype: pass "NA" to skip that level
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
cr_path <- args[1]; md_path <- args[2]
lineage <- args[3]; subtype <- args[4]; sub_subtype <- args[5]; out_csv <- args[6]

cr <- as.data.table(readRDS(cr_path))
md <- readRDS(md_path)

sel <- md$lineage == lineage
if (subtype != "NA") sel <- sel & md$lineage_subtype == subtype
if (sub_subtype != "NA") sel <- sel & md$lineage_sub_subtype == sub_subtype
lin_ids <- md$depmap_id[sel]
lin_set <- setdiff(unique(lin_ids), NA)
cat("lineage lines selected:", length(lin_set), "\n")

cat("CRISPR rows:", nrow(cr), "\n")

agg <- cr[, .(dep = mean(dependency)), by = .(gene_name, depmap_id)]
cat("agg rows:", nrow(agg), "\n")

wide <- dcast(agg, depmap_id ~ gene_name, value.var = "dep")
mat <- as.matrix(wide[, -1]); rownames(mat) <- wide$depmap_id
cat("wide matrix:", dim(mat), "\n")

is_lin <- rownames(mat) %in% lin_set
x <- mat[is_lin, ]
y <- mat[!is_lin, ]

nx <- colSums(!is.na(x)); ny <- colSums(!is.na(y))
mx <- colMeans(x, na.rm = TRUE); my <- colMeans(y, na.rm = TRUE)
vx <- colSums((x - rep(mx, each = nrow(x)))^2, na.rm = TRUE) / pmax(nx - 1, 1)
vy <- colSums((y - rep(my, each = nrow(y)))^2, na.rm = TRUE) / pmax(ny - 1, 1)

keep <- nx >= 8 & ny >= 100
se <- sqrt(vx[keep]/nx[keep] + vy[keep]/ny[keep])
tstat <- (mx[keep] - my[keep]) / se
dfree <- (vx[keep]/nx[keep] + vy[keep]/ny[keep])^2 /
      ((vx[keep]/nx[keep])^2/(nx[keep]-1) + (vy[keep]/ny[keep])^2/(ny[keep]-1))
pval <- 2 * pt(-abs(tstat), df = dfree)

tt <- data.frame(gene = names(mx[keep]), n_lineage = nx[keep], mean_lineage = mx[keep],
                 mean_rest = my[keep], delta = mx[keep] - my[keep],
                 t = tstat, p = pval, stringsAsFactors = FALSE)
tt$q <- p.adjust(tt$p, method = "BH")

frac_strong <- colMeans(mat < -0.5, na.rm = TRUE)
tt$frac_strong_all <- frac_strong[match(tt$gene, names(frac_strong))]

# Selective dependency: mean Chronos < -0.5, delta < -0.25, BH q < 0.05,
# not a common essential (strong effect in < 50% of all lines)
tt$selective <- with(tt, mean_lineage < -0.5 & delta < -0.25 & q < 0.05 & frac_strong_all < 0.5)
tt <- tt[order(tt$mean_lineage), ]
write.csv(tt, out_csv, row.names = FALSE)
cat("Selective dependencies:", sum(tt$selective), "of", nrow(tt), "\n")
