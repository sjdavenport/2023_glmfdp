# install.packages("BiocManager")
# BiocManager::install("GEOquery")
library("GEOquery")
library("stringr")
library("dplyr")
library("forcats")

# load series and platform data from GEO
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
gset <- getGEO("GSE42057", GSEMatrix = TRUE, getGPL = FALSE)
gset <- gset[[1]]

Y <- exprs(gset) # already log2 transformed

pheno <- phenoData(gset)
pdata <- pheno@data


# subject/control status
group_chr <- pdata[["source_name_ch1"]]
table(group_chr)
# Control subject PBMC    COPD subject PBMC 
# 42                   94 

group <- stringr::str_replace(pheno@data[["source_name_ch1"]], " subject PBMC", "")
covar_data <- pheno@data %>% 
    select(ends_with(":ch1")) %>% 
    rename_with(~gsub(":ch1", "", .x, fixed = TRUE))
str(covar_data)

# drop useless variables
to_drop <- c("tissue", "pctemph_slicer", "pctgastrap_slicer")
covar_data <- covar_data %>%
    select(!all_of(to_drop)) 

# export to txt files
X <- cbind(group = group, 
           covar_data)
str(X)
str(rownames(X))

write.table(X, file = "GSE42057_X.csv", sep = ",", quote = FALSE)
write.table(Y, file = "GSE42057_Y.csv", sep = ",")