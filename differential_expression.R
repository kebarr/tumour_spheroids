library(tidyr)
library(tibble)
library("DEP")
library("dplyr")
library(SummarizedExperiment)
library(ggplot2)
library(clusterProfiler)
library(collections)
library(org.Hs.eg.db)

setwd("/Users/user/Documents/proteomics")

AYA_GO<-read.csv("proteomics_gene_names.csv") # progenesis file, with awk used to get the gene names, and paste to add these as an extra column

label <- c('u87_Go100.1.1', 'u87_Go100.2.1', 'u87_Go100.3.1', 'u87_naive1.1', 'u87_naive2.1', 'u87_naive3.1')
condition <- c('u87_Go100.1', 'u87_Go100.1', 'u87_Go100.1', 'u87_naive.1', 'u87_naive.1', 'u87_naive.1')
replicate <- c('1', '2', '3','1', '2', '3')


exp_design <- data.frame(label, condition, replicate)
AYA_GO_unique <- make_unique(AYA_GO, "Gene.name", "X.Accession", delim = ";")
u87_columns <- grep("u87", colnames(AYA_GO_unique))

# DEP functions throw errors about data format, these work for our data
make_se <- function(proteins_unique, columns, expdesign) {
        raw <- proteins_unique[, columns]
        raw[raw == 0] <- NA
        raw <- log2(raw)
        expdesign <- mutate(expdesign, condition = make.names(condition)) %>% unite(ID, condition, replicate, remove = FALSE)
        rownames(expdesign) <- expdesign$ID
        matched <- match(make.names(expdesign$label), make.names(colnames(raw)))
        colnames(raw)[matched] <- expdesign$ID
        raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
        row_data <- proteins_unique[, -columns]
        rownames(row_data) <- row_data$name
        se <- SummarizedExperiment(assays = as.matrix(raw),colData = expdesign, rowData = row_data)
    return(se)
}


filter_missval <- function(se, thr=0, max_repl=3){
        thr <- 0
        max_repl <-3
        bin_data <- assay(se)
        idx <- is.na(assay(se))
        bin_data[!idx] <- 1
        bin_data[idx] <- 0
        keep <- bin_data %>%
            data.frame() %>%
            rownames_to_column() %>%
            gather(ID, value, -rowname) %>%
            left_join(., data.frame(colData(se)), by = "ID") %>%
            group_by(rowname, condition) %>%
            summarize(miss_val = n() - sum(value)) %>%
            filter(miss_val <= 0) %>%
            spread(condition, miss_val)
        se_fltrd <- se[keep$rowname, ]
    return(se_fltrd)
}


# steps as outline in Introduction to DEP (https://bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html#generate-a-data.frame-from-the-resulting-summarizedexperiment-object)
data_se <-make_se(AYA_GO_unique, u87_columns[7:12], exp_design)
data_filt <-filter_missval(data_se)
data_norm <- normalize_vsn(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_diff <- test_diff(data_imp, type = "control", control = "u87_naive.1")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

plot_volcano(dep, contrast = "u87_Go100.1_vs_u87_naive.1", label_size = 4, add_names = TRUE)

data_results <- get_results(dep)

sig <- data_results %>% filter(significant) # 37 significant, ratio column is average log 2 fold change
for_output <- c('name', 'ID', 'u87_Go100.1_vs_u87_naive.1_p.val', 'u87_Go100.1_vs_u87_naive.1_p.adj', 'u87_Go100.1_vs_u87_naive.1_ratio')
sig_output<-sig[, names(sig) %in% for_output]
write.csv(sig_output, "significant_dep_proteins.csv")


# use the significantly differentially expressed proteins for GO analysis

prots <- sig$ID

# get IDs in format required for GO function

genes <-mapIds(org.Hs.eg.db, as.character(prots), "ENTREZID", "UNIPROT")

gene_id_dict = Dict$new(genes)
 #populate list of ensemble names to pass into enrichGO
 l_entrez = c()
 for (i in prots){
    l_entrez = c(l_entrez, gene_id_dict$get(i))
    }

# perform GO analysis and plot results:


go_cc <- enrichGO(gene          = l_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

go_mf <- enrichGO(gene          = l_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

go_bp <- enrichGO(gene          = l_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

p_cc <- dotplot(go_cc, orderBy="GeneRatio")
p_mf <- dotplot(go_mf, orderBy="GeneRatio")
p_bp <- dotplot(go_bp, orderBy="GeneRatio")
