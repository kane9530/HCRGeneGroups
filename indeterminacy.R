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
processcsv_hes_cdh <- function(my_df, gene_names, my_sample_hes_cdh){
  my_df$sample_ID = my_sample_hes_cdh
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

#navigating to the right cdh6 hes6 folder


dir_hes_cdh <- list.dirs(path = "./dataCsv/hes6_cdh6")[-1]

my_hes_cdh_csv<- dir_hes_cdh %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_hes_cdh <- c()
for (dir_hes_cdh in my_hes_cdh_csv){
  my_ID <- rm_between(dir_hes_cdh, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_hes_cdh<- c(my_sample_hes_cdh, my_ID)
}

my_sample_hes_cdh <- as.integer(my_sample_hes_cdh)

gene_names = c("sox2", "tbxta","cdh6","hes6")

my_csv_files_cdh_hes <- lapply(my_hes_cdh_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

#merge sox2 and tbxta
my_csv_files_proc_hes_cdh_a <- vector("list", length = length(my_csv_files_cdh_hes))
my_csv_files_proc_merged_hes_cdh_a <- vector("list", length = length(my_csv_files_cdh_hes))

for (i in 1:length(my_csv_files_cdh_hes)){
  for (j in 1:length(gene_names)){
    my_csv_processed_cdh_hes_a <- processcsv_hes_cdh(my_df = my_csv_files_cdh_hes[[i]][[j]], 
                                        gene_names = gene_names[j],
                                        my_sample_hes_cdh = my_sample_hes_cdh[i])
    my_csv_files_proc_hes_cdh_a[[i]][[j]] <- my_csv_processed_cdh_hes_a
  }
  my_csv_files_proc_merged_hes_cdh_a[[i]] <- merge(my_csv_files_proc_hes_cdh_a[[i]][[1]], my_csv_files_proc_hes_cdh_a[[i]][[2]], by = c("cell_ID", "sample_ID"))
  
}
# merge cdh5 and hes6
my_csv_files_proc_hes_cdh_b <- vector("list", length = length(my_csv_files_cdh_hes))
my_csv_files_proc_merged_hes_cdh_b <- vector("list", length = length(my_csv_files_cdh_hes))

for (i in 1:length(my_csv_files_cdh_hes)){
  for (j in 1:length(gene_names)){
    my_csv_processed_cdh_hes_b <- processcsv_hes_cdh(my_df = my_csv_files_cdh_hes[[i]][[j]], 
                                                     gene_names = gene_names[j+2],
                                                     my_sample_hes_cdh = my_sample_hes_cdh[i])
    my_csv_files_proc_hes_cdh_b[[i]][[j]] <- my_csv_processed_cdh_hes_b
  }
  my_csv_files_proc_merged_hes_cdh_b[[i]] <- merge(my_csv_files_proc_hes_cdh_b[[i]][[1]], my_csv_files_proc_hes_cdh_b[[i]][[2]], by = c("cell_ID", "sample_ID"))
}


#my_csv_files_proc_merged_hes_cdh <- vector("list", length = length(my_csv_files_cdh_hes))


df_combined_a <- my_csv_files_proc_merged_hes_cdh_a %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])


df_combined_b <- my_csv_files_proc_merged_hes_cdh_b %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) 
  

df_combined <- merge(df_combined_a, df_combined_b, by = c("cell_ID", "sample_ID")) %>% 
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 8:9) %>%
  ungroup()%>%
  mutate(lineage_index = .[[10]]-.[[11]])

#plotting

p_nmps_only <- ggplot(df_combined, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth = 0.2, color="black")+
  scale_x_continuous(name = "NMindex",limits=c(-1, 1))+
  theme_minimal()+
  scale_fill_manual(values=plasma(n=length(my_csv_files_cdh_hes)))+
  theme(text = element_text(size = 15, family = "sans"))+
  ggtitle("NMP distribution")+
  #theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))


line_nmps_only <- ggplot(df_combined, aes(x=nm_index, y=))+
  geom_density(colour = "red")+
  theme_minimal()+
  ggtitle("NMP distribution")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size = 15, family = "sans"))

# plotting the lineage marker index
p_lineage_only <- ggplot(df_combined, aes(x=lineage_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth = 0.2, color="black")+
  scale_x_continuous(name = "Lineage index",limits=c(-1, 1))+
  theme_minimal()+
  scale_fill_manual(values=plasma(n=length(my_csv_files_cdh_hes)))+
  theme(text = element_text(size = 15, family = "sans"))+
  ggtitle("Lineage marker distribution")+
  #theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

p_lineage_only

line_lineage <- ggplot(df_combined, aes(x=lineage_index, y=))+
  geom_density(colour = "blue")+
  theme_minimal()+
  ggtitle("Lineage marker distribution")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size = 15, family = "sans"))

line_lineage

#plotting both nm index and lineage

p_nmp_lineage <- ggplot(df_combined, aes(x=nm_index, y=lineage_index))+
  geom_point(color = "grey")+
  geom_smooth(method = "lm", color="navy", fill="lightblue")+
  theme_minimal()+
  labs(x = "NM index", y = "Lineage index", title = "Indeterminacy of NMPs")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size = 15, family = "sans"))+
  stat_cor(method = "pearson", label.x = -1, label.y = 0.6)
p_nmp_lineage

#calculate pearsons coefficient
# use options(scipen=0) to force different scientific notation
PCC <- cor.test(df_combined$nm_index, df_combined$lineage_index, method = "pearson", exact = T, conf.level = 0.95)

PCC
PCC$estimate
PCC$conf.int[c(1,2)]
PCC$p.value

