# Genome-wide per-gene Cox models on TCGA tumors (config-driven)
# Usage: Rscript survival_cox.R <expr_csv.gz> <survival_csv> <clinical_csv> <probemap> <out_csv>
# Expression continuous (z-scored log2 TPM), covariates: age, late stage (>= III)
suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = TRUE)
expr_path <- args[1]; surv_path <- args[2]; clin_path <- args[3]; pm_path <- args[4]
out_csv <- args[5]

expr <- read.csv(gzfile(expr_path), row.names = 1, check.names = FALSE)
mat <- as.matrix(expr)

pm <- read.delim(pm_path, stringsAsFactors = FALSE)
map <- pm$gene[match(rownames(mat), pm$id)]
ok <- !is.na(map); mat <- mat[ok, ]; map <- map[ok]
mn <- rowMeans(mat); ord <- order(map, -mn); mat <- mat[ord, ]; map <- map[ord]
d <- !duplicated(map); mat <- mat[d, ]; rownames(mat) <- map[d]

tum <- mat[, grepl("-0[16]$", colnames(mat))]
patients <- sub("^TCGA-(..)-(....).*$", "\\1-\\2", colnames(tum))
colnames(tum) <- patients

surv <- read.csv(surv_path, row.names = 1, check.names = FALSE)
surv <- as.data.frame(t(surv))
surv$patient <- sub("^TCGA-(..)-(....).*$", "\\1-\\2", rownames(surv))

clin <- read.csv(clin_path, row.names = 1, check.names = FALSE)
clin <- as.data.frame(t(clin))
clin$patient <- sub("^TCGA-(..)-(....).*$", "\\1-\\2", rownames(clin))

age <- suppressWarnings(as.numeric(clin$age_at_index.demographic))
stage_raw <- clin$ajcc_pathologic_stage.diagnoses
names(age) <- clin$patient; names(stage_raw) <- clin$patient

df <- data.frame(patient = surv$patient,
                 OS.time = as.numeric(surv[["OS.time"]]),
                 OS = as.numeric(surv[["OS"]]),
                 age = age[surv$patient],
                 stringsAsFactors = FALSE)

stage_num <- suppressWarnings(as.numeric(stage_raw[df$patient]))
df$late_stage <- ifelse(is.na(stage_num), NA, ifelse(stage_num >= 3, 1, 0))

df <- df[!is.na(df$OS.time) & !is.na(df$OS) & df$OS.time > 0, ]
cat("Patients with survival data:", nrow(df), "\n")

common <- intersect(df$patient, colnames(tum))
cat("With expression:", length(common), "\n")
tum <- tum[, common]
df <- df[match(common, df$patient), ]

expr_ok <- rowSums(tum > 1) >= 0.1 * ncol(tum)
tum <- tum[expr_ok, ]
cat("Expressed genes for survival testing:", sum(expr_ok), "\n")

z <- t(scale(t(tum)))

res <- do.call(rbind, lapply(rownames(z), function(g) {
  d <- df
  d$exprz <- z[g, d$patient]
  ok <- complete.cases(d[, c("OS.time", "OS", "exprz", "age", "late_stage")])
  d <- d[ok, ]
  if (nrow(d) < 80) return(NULL)
  f <- coxph(Surv(OS.time, OS) ~ exprz + age + late_stage, data = d)
  s <- summary(f)
  data.frame(gene = g, n = nrow(d), HR = s$coefficients["exprz", "exp(coef)"],
             p = s$coefficients["exprz", "Pr(>|z|)"],
             coef = s$coefficients["exprz", "coef"], stringsAsFactors = FALSE)
}))
res$q <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
write.csv(res, out_csv, row.names = FALSE)
cat("Significant (FDR<0.05):", sum(res$q < 0.05), "\n")
