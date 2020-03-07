#loading libraries

library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(qdapRegex)
library(viridis)
library(ggpubr)
library(extrafont)

#setting up processing functions, ripped from Kane's previous code
processcsv_hes_cdh_1 <- function(my_df, gene_names, my_sample_hes_cdh_1){
  my_df$sample_ID = my_sample_hes_cdh_1
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

#navigating to the right cdh6 hes6 folder


dir_hes_cdh <- list.dirs(path = "./dataCsv/sox2_tbxta_hes6_cdh6")[-1]

my_hes_cdh_csv<- dir_hes_cdh %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_no_18 <- c()
for (dir18ss in my_dir18ss_csv){
  my_ID_18 <- rm_between(dir18ss, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_no_18 <- c(my_sample_no_18, my_ID_18)
}

my_sample_no_18 <- as.integer(my_sample_no_18)

gene_names = c("sox2", "tbxta")

my_csv_files_18 <- lapply(my_dir18ss_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_18 <- vector("list", length = length(my_csv_files_18))
my_csv_files_proc_merged_18 <- vector("list", length = length(my_csv_files_18))

for (i in 1:length(my_csv_files_18)){
  for (j in 1:length(gene_names)){
    my_csv_processed_18 <- processcsv18(my_df = my_csv_files_18[[i]][[j]], 
                                        gene_names = gene_names[j],
                                        my_sample_no_18 = my_sample_no_18[i])
    my_csv_files_proc_18[[i]][[j]] <- my_csv_processed_18
  }
  my_csv_files_proc_merged_18[[i]] <- merge(my_csv_files_proc_18[[i]][[1]], my_csv_files_proc_18[[i]][[2]], by = c("cell_ID", "sample_ID"))
  #my_csv_files_proc_merged_18[[i]] <- bind_rows(list(my_csv_files_proc_18[i][1],my_csv_files_proc_18[i][2],my_csv_files_proc_18[i][3]))
}


df_combined_18 <- my_csv_files_proc_merged_18 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#plotting

