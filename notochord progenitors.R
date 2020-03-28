#Importing all libraries I could possibly want to use

library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(qdapRegex)
library(viridis)
library(ggpubr)
library(extrafont)
library(ggridges)
library(RColorBrewer)


#####processing for notochord progenitors
#####Where I say mp I really mean notochord progenitor, it's just that I didn't know they were NPs
#####Will tidy up the code when I have submitted my project report haha

processcsvmp <- function(my_df, gene_names, my_sample_mp){
  my_df$sample_ID = my_sample_mp
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}


##### Importing MP data

dirmp <- list.dirs(path = "./dataCsv/notochord_progenitors")[-1]

my_dir_mp_csv<- dirmp %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_mp <- c()
for (dirmp in my_dir_mp_csv){
  my_ID_mp <- rm_between(dirmp, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_mp <- c(my_sample_mp, my_ID_mp)
}

my_sample_mp <- as.integer(my_sample_mp)

gene_names = c("sox2", "tbxta")

my_csv_files_mp <- lapply(my_dir_mp_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_mp<- vector("list", length = length(my_csv_files_mp))
my_csv_files_proc_merged_mp <- vector("list", length = length(my_csv_files_mp))

for (i in 1:length(my_csv_files_mp)){
  for (j in 1:length(gene_names)){
    my_csv_processed_mp <- processcsvmp(my_df = my_csv_files_mp[[i]][[j]], 
                                        gene_names = gene_names[j],
                                        my_sample_mp = my_sample_mp[i])
    my_csv_files_proc_mp[[i]][[j]] <- my_csv_processed_mp
  }
  my_csv_files_proc_merged_mp[[i]] <- merge(my_csv_files_proc_mp[[i]][[1]], my_csv_files_proc_mp[[i]][[2]], by = c("cell_ID", "sample_ID"))
  
}

##### Turning data into dataframe

df_combined_mp <- my_csv_files_proc_merged_mp %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

###### Plotting into graphs

p_mp <- ggplot(df_combined_mp, aes(x=nm_index))+ 
  geom_histogram()+
  labs(title = "Notochord progenitors")
p_mp

##### What happens if I normalise according to NMP population?

max_sox2 <- 903
min_sox2 <- 6

max_tbxta <- 1118.11
min_tbxta <- 6

df_combined_mesodermal_progenitors <- my_csv_files_proc_merged_mp %>% 
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(sox2_normalise = ~((.-min_sox2)/(max_sox2-min_sox2))), .vars = 3) %>%
  mutate_at(.funs = list(tbxta_normalise = ~((.-min_tbxta)/(max_tbxta-min_tbxta))), .vars = 4)%>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

### Plotting some NMP normalised mesodermal progenitors

p_mp_norm_nmp <- ggplot(df_combined_mesodermal_progenitors, aes(x=nm_index))+
  geom_histogram()+
  labs(title = "Notochord progenitors normalised to NMPS")

p_mp_norm_nmp

line_mp_norm_nmp <- ggplot(df_combined_mesodermal_progenitors, aes(x=nm_index))+
  geom_density()+
  labs(title = "Notochord progenitors normalised to NMPS")+
  theme_minimal()

line_mp_norm_nmp

##### Importing the neural tube data here and normalising


## processing function for neural tube cells
processcsvneural <- function(my_df, gene_names, my_sample_neural){
  my_df$sample_ID = my_sample_neural
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

##### The neural cells #####

## navigate to the neural tube directory

dirneural <- list.dirs(path = "./dataCsv/neural_tube")[-1]

my_dir_neural_csv<- dirneural %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_neural <- c()
for (dirneural in my_dir_neural_csv){
  my_ID_neural <- rm_between(dirneural, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_neural <- c(my_sample_neural, my_ID_neural)
}

my_sample_neural <- as.integer(my_sample_neural)

gene_names = c("sox2", "tbxta")

my_csv_files_neural <- lapply(my_dir_neural_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_neural <- vector("list", length = length(my_csv_files_neural))
my_csv_files_proc_merged_neural <- vector("list", length = length(my_csv_files_neural))

for (i in 1:length(my_csv_files_neural)){
  for (j in 1:length(gene_names)){
    my_csv_processed_neural <- processcsvneural(my_df = my_csv_files_neural[[i]][[j]], 
                                                gene_names = gene_names[j],
                                                my_sample_neural = my_sample_neural[i])
    my_csv_files_proc_neural[[i]][[j]] <- my_csv_processed_neural
  }
  my_csv_files_proc_merged_neural[[i]] <- merge(my_csv_files_proc_neural[[i]][[1]], my_csv_files_proc_neural[[i]][[2]], by = c("cell_ID", "sample_ID"))
  
}

df_combined_neural <- my_csv_files_proc_merged_neural %>% 
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID)  %>%
  mutate_at(.funs = list(sox2_normalise = ~((.-min_sox2)/(max_sox2-min_sox2))), .vars = 3) %>%
  mutate_at(.funs = list(tbxta_normalise = ~((.-min_tbxta)/(max_tbxta-min_tbxta))), .vars = 4)%>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

###According to my quantification sheet, only samples 5,6,7,8,10 and 11 should be used for the NT

df_combined_neural_filtered <- df_combined_neural %>% 
  filter(sample_ID == 5|
         sample_ID == 6|
         sample_ID == 7|
         sample_ID == 8|
         sample_ID == 10|
         sample_ID == 11)

##Quick check

plot_neural_filtered <- ggplot(df_combined_neural_filtered, aes(x=nm_index))+
  geom_histogram()
plot_neural_filtered

plot_neural_unfiltered <- ggplot(df_combined_neural, aes(x=nm_index))+
  geom_histogram()
plot_neural_unfiltered


## Combining dataframes


df_combined_het <- bind_rows(df_combined_mesodermal_progenitors, df_combined_neural_filtered, .id = "origin")

df_combined_het$origin <- as.factor(df_combined_het$origin)

levels(df_combined_het$origin) <- c("NP", "NT")

#Plotting

line_both <- ggplot(df_combined_het, aes(x=nm_index, colour=factor(origin)))+
  geom_density(size=2)+
  scale_x_continuous( limits=c(-0.75, 0.75))+
  guides(colour=guide_legend(title = "Type of cell"))+
  theme_minimal()+
  labs(x="NM Index", y="Density", title = "NM index in Notochord Progenitors and Neural Tube cells")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = viridis(n=2, begin = 0.05, end = 0.95))

line_both

##### Importing NMP data so that I can plot it all into a graph


##Processing function for 18ss

processcsv18 <- function(my_df, gene_names, my_sample_no_18){
  my_df$sample_ID = my_sample_no_18
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

#####The normal 18ss NMPS ######

##Navigating to the right 18ss directory

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

#Converting 18ss into dataframe
df_combined_18 <- my_csv_files_proc_merged_18 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

##### Making an unholy behemoth of a dataframe, combining 18ss with neural and mesodermal data

df_all_3 <- bind_rows(df_combined_mesodermal_progenitors, df_combined_neural_filtered, df_combined_18, .id = "origin" )

df_all_3$origin <- as.factor(df_all_3$origin)

levels(df_all_3$origin) <- c("NP", "NT", "NMP")

### Plotting all 3 populations of cells into a graph

line_all_3 <- ggplot(df_all_3, aes(x=nm_index, fill=factor(origin)))+
  geom_histogram(position = "identity", alpha = 0.5, colour="grey2")+
  theme_minimal()+
  labs(x="NM index", y="Density")+
  ggtitle("Comparisons of heterogeneity across 3 populations of cells")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend("Type of cell"))+
  scale_x_continuous(limits = c(-1, 1))+
  #scale_fill_manual(values = c("orange", "yellow", "red"))
  scale_fill_manual(values = viridis(n=3))+
  theme(text = element_text(size = 15))
line_all_3

plot_all_3 <- ggplot(df_all_3, aes(x=nm_index, fill=factor(origin)))+
  geom_density()
plot_all_3
### Trying out a new plotting function, ggridges

df_all_3$origin <- factor(df_all_3$origin, levels = c("NMP", "NT", "NP"))

ridge_all_3 <- ggplot(df_all_3, aes(x=nm_index, y=origin, fill = factor(origin)))+
  geom_density_ridges(scale=5, show.legend = T, alpha = 0.5)+
  theme_minimal()+
  scale_fill_manual(values = viridis(n=3, direction = -1, begin = 0.05, end = 0.95))+
  scale_x_continuous(limits = c(-1, 1))+
  labs(x="NM index", y="Origin of cell")+
  ggtitle("The heterogeneity of different cell populations")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend("Type of Cell"))+
  theme(text = element_text(size = 15))

ridge_all_3

