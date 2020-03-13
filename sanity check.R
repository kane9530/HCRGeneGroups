#sanity check
#loading libraries

library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(qdapRegex)
library(viridis)
library(ggpubr)
library(extrafont)

# 18ss
processcsv18 <- function(my_df, gene_names, my_sample_no_18){
  my_df$sample_ID = my_sample_no_18
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

dir18ss <- list.dirs(path = "./dataCsv/18ss")[-1]

my_dir18ss_csv<- dir18ss %>% 
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
  
}


df_combined_18 <- my_csv_files_proc_merged_18 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#Looking at specifically hes6 and cdh6 datasets of nmps

df_combine_18_hes6_cdh6_only <- df_combined_18 %>% 
  filter(sample_ID == 10|sample_ID == 11|sample_ID == 12|sample_ID == 13|sample_ID == 14)

#View(df_combined_18)

p_18 <- ggplot(df_combine_18_hes6_cdh6_only, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth=0.2, color = 'black')+
  scale_x_continuous( limits=c(-1, 1))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=length(my_csv_files_18)))+
  theme_minimal() + 
  labs(x="NM index", y="Number of Cells")+
  ggtitle("NMPs of hes6 and cdh6 samples") + 
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))
p_18

line_18 <- ggplot(df_combine_18_hes6_cdh6_only, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+ 
  ggtitle("NMPs of hes6 and cdh6 samples")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")
line_18

#####Time to import the hes6 and cdh6 data here#####


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


df_combined$sample_ID <- as.factor(df_combined$sample_ID)
# Plots of NMPs that co-express Cdh6 and Hes6


p_nmps_co_only <- ggplot(df_combined, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth = 0.2, color="black")+
  scale_x_continuous(name = "NMindex",limits=c(-1, 1))+
  theme_minimal()+
  scale_fill_manual(values=plasma(n=length(my_csv_files_cdh_hes)))+
  theme(text = element_text(size = 15, family = "sans"))+
  ggtitle("NMP distribution")+
  #theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))

p_nmps_co_only


line_nmps_co_only <- ggplot(df_combined, aes(x=nm_index, y=))+
  geom_density(colour = "red")+
  theme_minimal()+
  ggtitle("NMP distribution")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size = 15, family = "sans"))

line_nmps_co_only

##Time to plot the number of NMPs vs the number of NMPS coexpressing Cdh6 and Hes6
#Will need to merge dataframes

merged_df_nmps <- df_combine_18_hes6_cdh6_only %>% 
  bind_rows(df_combined, .id = "id")

merged_df_nmps$id <- gsub('1','NMPs', merged_df_nmps$id)

merged_df_nmps$id <- gsub('2','Indeterminates',merged_df_nmps$id)

merged_df_nmps$id <- as.factor(merged_df_nmps$id)
##Plotting time

#The following graph doesn't mean much, it has some suggestions of a bell shaped curve
#with a bimodal distribution super-imposed on top
nmps_and_co_nmps <- ggplot(merged_df_nmps, aes(x=nm_index, y=, colour = factor(id)))+
  geom_density(size=2)+
  guides(colour=guide_legend(title="Type of cell"))+
  theme_minimal()

nmps_and_co_nmps

numbers_of_nmps_vs_co_nmps <-  ggplot(merged_df_nmps, aes(x=id, y=, fill = factor(id)))+
  geom_bar()+
  ggtitle("Number of NMPs vs number of All-4's")+
  theme_minimal()+
  guides(fill=guide_legend(title="Type of cell"))
numbers_of_nmps_vs_co_nmps

##Kane asked for a breakdown of what constitutes the indeterminates.

further_breakdown

## Find the maximum cdh6 and maximum hes6's 
# Find maximum of sox2_normalise and tbxta_normalise
# Ask Kane for an appropriate threshold
# List the Cell ID and Sample ID
#filter


