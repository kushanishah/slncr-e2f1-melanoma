library(tidyverse)

# variables
exp = read_tsv("MMyeloma1nam.csv")
n = ncol(exp) - 2
pct = round(n * 0.33)
lncrna = "ENSG00000215417_10"
difference = 0.5
###########


exp_ma = exp[, 3:ncol(exp)] %>% as.matrix()
# keep expressed genes and with some variation. Change the numebrs as you wish
min_sum_exp = 5000 # change this as needed
min_sd_exp = 2 # change this as needed
exp_keep = rowSums(exp_ma)>min_sum_exp & apply(exp_ma, 1, sd)>min_sd_exp
sum(exp_keep) # this is the number of total genes. Anything above 15000 could be very slow and not end ever.

lncrna_exp = exp %>% filter(Description==lncrna) %>% 
    select(-Name) %>% 
    gather(Description, value) %>% 
    arrange(value)

# get top and bottom part
low = lncrna_exp %>% head(pct) %>% pull(Description)
top = lncrna_exp %>% tail(pct) %>% pull(Description)

low_ma = exp[exp_keep,c("Name", low)] %>% 
    column_to_rownames("Name") %>% 
    as.matrix()
top_ma = exp[exp_keep,c("Name", top)] %>% 
    column_to_rownames("Name") %>% 
    as.matrix()

# correlation function
fast_cor = function(ma){
    maxn = ifelse(nrow(ma)<70000, nrow(ma), 70000) # change 15000 to the number of genes you want to limit the correlation to
    message("Using: ", maxn, " genes.")
    cor_ma = cor(t(ma[1:maxn,]), method = "spearman")
    cor_ma
}

low_cor = fast_cor(low_ma)
top_cor = fast_cor(top_ma)

d = low_cor - top_cor
keep_d = abs(d) > difference
k = arrayInd(which(keep_d), dim(d))
valid = data.frame(gene_1=rownames(d)[k[,1]],
                   gene_2=colnames(d)[k[,2]],
                   cor_low=low_cor[keep_d],
                   cor_top=top_cor[keep_d],
                   difference=d[keep_d])
valid = valid %>% distinct() %>% 
    left_join(exp[,1:2], by = c("gene_1" = "Name")) %>% 
    left_join(exp[,1:2], by = c("gene_2" = "Name"))
write_csv(valid, "1.csv") # change name  file


## select PARTNERS
partners <- readxl::read_xlsx("alllmir17Hgpartnerensgand.xlsx", sheet = 1, col_names = FALSE)[,1, drop = TRUE]
valid <- valid %>%
    filter(Description.x %in% partners | Description.y %in% partners)
write_csv(valid, "2.csv")

## select TARGETS
targets <- readxl::read_xlsx("readyforwertp_targetsofmir17Hg.xlsx", sheet = 1, col_names = FALSE)[,1, drop = TRUE]
valid <- valid %>%
    filter(Description.x %in% targets | Description.y %in% targets)
write_csv(valid, "3.csv")
