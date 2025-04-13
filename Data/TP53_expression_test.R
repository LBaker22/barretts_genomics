
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery") 

install.packages("ggplot2")



library(GEOquery)
library(limma)
library(ggplot2)

# 1. Download GEO data
gse <- getGEO("GSE26886", GSEMatrix = TRUE)
eset <- gse[[1]]

# 2. Extract phenotype data
pheno <- pData(eset)
table(pheno$characteristics_ch1.1)  # Check sample groups

# 3. Create a simplified grouping
pheno$group <- ifelse(grepl("Barrett", pheno$characteristics_ch1), "Barrett",
                      ifelse(grepl("normal", pheno$characteristics_ch1), "Normal", NA))

# 4. Filter for Barrett and Normal samples
keep <- which(!is.na(pheno$group))
eset_sub <- eset[, keep]
pheno_sub <- pheno[keep, ]

# 5. Extract TP53 expression
exprs_mat <- exprs(eset_sub)
tp53_probe <- grep("TP53", fData(eset_sub)$Gene.symbol, value = FALSE)
tp53_expr <- exprs_mat[tp53_probe, ]

# 6. Create a data frame for plotting
df <- data.frame(
  Expression = as.numeric(tp53_expr),
  Group = pheno_sub$group
)

# 7. Plot and test
ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  labs(title = "TP53 Expression in Barrett's vs Normal Oesophagus",
       y = "TP53 Expression (log2 intensity)") +
  theme_minimal()

# 8. Statistical test
t_test <- t.test(Expression ~ Group, data = df)
print(t_test)
