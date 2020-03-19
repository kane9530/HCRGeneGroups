##### Trying to justify the heterogeneity of NMPs vs noto and neural tube ######

library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(qdapRegex)
library(viridis)
library(ggpubr)
library(extrafont)


##Processing function for 18ss

processcsv18 <- function(my_df, gene_names, my_sample_no_18){
  my_df$sample_ID = my_sample_no_18
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

## processing function for neural tube cells
processcsvneural <- function(my_df, gene_names, my_sample_neural){
  my_df$sample_ID = my_sample_neural
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}


## processing function for notochord cells

processcsvnoto <- function(my_df, gene_names, my_sample_noto){
  my_df$sample_ID = my_sample_noto
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

#####The notochord cells#####

##navigating to the notochord data

dir_noto <- list.dirs(path = "./dataCsv/notochord")[-1]

my_dir_noto_csv<- dir_noto %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_noto <- c()
for (dir_noto in my_dir_noto_csv){
  my_ID_noto <- rm_between(dir_noto, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_noto <- c(my_sample_noto, my_ID_noto)
}

my_sample_noto <- as.integer(my_sample_noto)

gene_names = c("sox2", "tbxta")

my_csv_files_noto <- lapply(my_dir_noto_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_noto <- vector("list", length = length(my_csv_files_noto))
my_csv_files_proc_merged_noto <- vector("list", length = length(my_csv_files_noto))


for (i in 1:length(my_csv_files_noto)){
  for (j in 1:length(gene_names)){
    my_csv_processed_noto <- processcsvnoto(my_df = my_csv_files_noto[[i]][[j]], 
                                            gene_names = gene_names[j],
                                            my_sample_noto = my_sample_noto[i])
    my_csv_files_proc_noto[[i]][[j]] <- my_csv_processed_noto
  }
  my_csv_files_proc_merged_noto[[i]] <- merge(my_csv_files_proc_noto[[i]][[1]], my_csv_files_proc_noto[[i]][[2]], by = c("cell_ID", "sample_ID"))
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

##### Processing all of the above data into a dataframe

##I will normalise in accordance with the NMPs

#Converting 18ss into dataframe
df_combined_18 <- my_csv_files_proc_merged_18 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#Find the minimum and maximum NMP sox2 and tbxta values to normalise notochord cells and neural tube cells
#to.

max_sox2 <- max(df_combined_18$sox2)
min_sox2 <- min(df_combined_18$sox2)

max_tbxta <- max(df_combined_18$tbxta)
min_tbxta <- min(df_combined_18$tbxta)

#Converting noto cells into dataframe and normalising with regards to the NMP population
df_combined_noto <- my_csv_files_proc_merged_noto %>% 
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(sox2_normalise = ~((.-min_sox2)/(max_sox2-min_sox2))), .vars = 3) %>%
  mutate_at(.funs = list(tbxta_normalise = ~((.-min_tbxta)/(max_tbxta-min_tbxta))), .vars = 4)%>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])
## notochord cells normalised to itself

df_combined_noto_self_norm <- my_csv_files_proc_merged_noto %>% 
  reduce(rbind) %>% 
  group_by(sample_ID) %>% 
  arrange(sample_ID, cell_ID) %>% 
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>% 
  ungroup() %>% 
  mutate(nm_index = .[[5]]-.[[6]])

#Converting neural cells into dataframe and normalising with regards to the NMP population
df_combined_neural <- my_csv_files_proc_merged_neural %>% 
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID)  %>%
  mutate_at(.funs = list(sox2_normalise = ~((.-min_sox2)/(max_sox2-min_sox2))), .vars = 3) %>%
  mutate_at(.funs = list(tbxta_normalise = ~((.-min_tbxta)/(max_tbxta-min_tbxta))), .vars = 4)%>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

##neural cells normalised to itself

df_combined_neural_self_norm <- my_csv_files_proc_merged_neural %>% 
  reduce(rbind) %>% 
  group_by(sample_ID) %>% 
  arrange(sample_ID, cell_ID) %>% 
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>% 
  ungroup() %>% 
  mutate(nm_index = .[[5]]-.[[6]])
  
#Merging all above dataframes together 

df_combined <- df_combined_18 %>% 
  bind_rows(df_combined_noto, df_combined_neural, .id = "origin")

#Renaming the origin column so that they reflect where the row comes from
df_combined$origin <- as.factor(df_combined$origin)

levels(df_combined$origin) <- c("nmps", "notochord","neural_tube")

#Check the dataframe looks normal
View(df_combined)

## Looking at the min and max values of noto and heural normalised

min_noto_sox2 <- min(df_combined_noto$sox2_normalise)
max_noto_sox2 <- max(df_combined_noto$sox2_normalise)

min_noto_tbxta <- min(df_combined_noto$tbxta_normalise)
max_noto_tbxta <- max(df_combined_noto$tbxta_normalise)

min_neural_sox2 <- min(df_combined_neural$sox2_normalise)
max_neural_sox2 <- max(df_combined_neural$sox2_normalise)

min_neural_tbxta <- min(df_combined_neural$tbxta_normalise)
max_neural_tbxta <- max(df_combined_neural$tbxta_normalise)

notochord_range_sox2 <- ggplot(df_combined_noto, aes(x=sox2_normalise))+
  geom_histogram(bindwidth=0.2)
notochord_range_sox2

notochord_range_tbxta <- ggplot(df_combined_noto, aes(x=tbxta_normalise))+
  geom_histogram(bindwidth=0.2)
notochord_range_tbxta

neural_range_tbxta_un_norm <- ggplot(df_combined_neural, aes(x=tbxta))+
  geom_histogram(binwidth = 0.5)
neural_range_tbxta_un_norm

sox2_un_norm <- ggplot(df_combined, aes(x=sox2, fill=factor(origin)))+
  geom_histogram(position= "identity", bindwidth=0.5, alpha = 0.5)
sox2_un_norm

tbxta_un_norm <- ggplot(df_combined, aes(x=tbxta, fill=factor(origin)))+
  geom_histogram(position= "identity", binwidth = 30, alpha = 0.5)
tbxta_un_norm

sox2_norm <- ggplot(df_combined, aes(x=sox2_normalise, fill=factor(origin)))+
  geom_histogram(position= "identity", bindwidth=0.5, alpha = 0.5)
sox2_norm

tbxta_norm <- ggplot(df_combined, aes(x=tbxta_normalise, fill=factor(origin)))+
  geom_histogram(position= "identity", bindwidth=0.5, alpha = 0.5)
tbxta_norm
#########sample 1 seems very off in tbxta values
###

neural_un_norm_max_sox2 <- max(df_combined_neural$sox2)
neural_un_norm_min_sox2 <- min(df_combined_neural$sox2)

nmp_un_norm_max_sox2 <- max(df_combined_18$sox2)
nmp_un_norm_min_sox2 <- min(df_combined_18$sox2)

neural_un_norm_max_tbxta <-max(df_combined_neural$tbxta)
neural_un_norm_min_tbxta <- min(df_combined_neural$tbxta)

nmp_un_norm_max_tbxta <- max(df_combined_18$tbxta)
nmp_un_norm_min_tbxta <- min(df_combined_18$tbxta)



#####Plotting NM index

# Plots of notochord cells

p_noto<- ggplot(df_combined_noto, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth=0.1, color = 'black')+
  scale_x_continuous( limits=c(-0.5, 0.5))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=length(my_csv_files_noto)))+
  theme_minimal() +
  labs(x="NM index", y="Number of Cells")+
  ggtitle("Notochords") +
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))

p_noto

line_noto <- ggplot(df_combined_noto, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+
  ggtitle("Notochords")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")
line_noto

##plotting notochords without sample 3


df_combined_noto_w_sample_1 <- df_combined_noto %>% 
  filter(sample_ID==1)

p_noto_w_sample_7<- ggplot(df_combined_noto_w_sample_7, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(bins = 3, color = 'black')+
  scale_x_continuous( limits=c(-0.5, 0.5))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=5))+
  theme_minimal() +
  labs(x="NM index", y="Number of Cells")+
  ggtitle("Notochords restricted to sample 7") +
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))

p_noto_w_sample_7

line_noto_w_sample_1 <- ggplot(df_combined_noto_w_sample_1, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+
  ggtitle("Notochords restricted to sample 4")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")

line_noto_w_sample_1

p_noto_s_1 <- ggplot(df_combined_noto_w_sample_1, aes(x=tbxta_normalise))+
  geom_histogram()
p_noto_s_1

df_combined_18_w_sample_2 <- df_combined_18 %>% 
  filter(sample_ID==2)

p_nmps_s_1 <- ggplot(df_combined_18_w_sample_2, aes(x=tbxta_normalise))+
  geom_histogram()
p_nmps_s_1


# Plots of neural tube cells

p_neural<- ggplot(df_combined_neural, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth=0.1, color = 'black')+
  scale_x_continuous( limits=c(-1.5, 1.5))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=length(my_csv_files_neural)))+
  theme_minimal() +
  labs(x="NM index", y="Number of Cells")+
  ggtitle("Neural Tube cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))

p_neural

line_neural <- ggplot(df_combined_neural, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+
  ggtitle("Neural Tube cells")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")

line_neural

# Plots of NMPs only

p_NMPs <- ggplot(df_combined_18,aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth=0.1, color = 'black')+
  scale_x_continuous( limits=c(-1.5, 1.5))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=length(my_csv_files_18)))+
  theme_minimal() +
  labs(x="NM index", y="Number of Cells")+
  ggtitle("NMPs") +
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))

p_NMPs

line_NMPs <- ggplot(df_combined_18, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+
  ggtitle("NMPs")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")

line_NMPs

#Plots of all 3 types of cells together

p_all_cells <- ggplot(df_combined, aes(x=nm_index, fill=factor(origin)))+
  geom_histogram(binwidth = 0.1, alpha = 0.20)+
  scale_x_continuous( limits=c(-1.5, 1.5))+
  theme(text=element_text(size = 30, family = "sans"))+
  theme_minimal() +
  labs(x="NM index", y="Number of Cells")+
  ggtitle("NM index in all cells") +
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Origin of cell"))

p_all_cells

line_all_cells <- ggplot(df_combined, aes(x=nm_index, fill=factor(origin)))+
  geom_density(size = 2, alpha = 0.5, show.legend = T, linetype = 0)+
  guides(fill=guide_legend(title="Origin of cell"))+
  theme_minimal()+
  scale_fill_manual(values = viridis(n=3))

line_all_cells

## looking at individual samples

sample_1 <- df_combined %>% 
  filter(sample_ID==1)
View(sample_1)

un_norm_sox2_s_1 <- ggplot(sample_1, aes(x=sox2, fill=factor(origin)))+
  geom_histogram(position= "identity", bindwidth=0.5, alpha = 0.5)
un_norm_sox2_s_1
####hmmmm, but the sample numbers don't actually line up across the files


sample_2 <- df_combined %>% 
  filter(sample_ID==2)
View(sample_2)

un_norm_sox2_s_2 <- ggplot(sample_2, aes(x=sox2, fill=factor(origin)))+
  geom_histogram(position= "identity", bindwidth=0.5, alpha = 0.5)
un_norm_sox2_s_2

##### plotting for self normalised notochord cells

p_tbxta_dist_self_norm_noto <- ggplot(df_combined_noto_self_norm, aes(x=tbxta_normalise))+
  geom_histogram()

p_tbxta_dist_self_norm_noto

p_tbxta_dist_not_norm_noto <- ggplot(df_combined_noto_self_norm, aes(x=tbxta))+
  geom_histogram()
p_tbxta_dist_not_norm_noto

p_tbxta_dist_not_norm_nmps <- ggplot(df_combined_18, aes(x=tbxta))+
  geom_histogram()
p_tbxta_dist_not_norm_nmps

##### plotting self normalised neural cells

p_nm_index_neural_self_norm <- ggplot(df_combined_neural_self_norm, aes(x=nm_index))+
  geom_histogram()

p_nm_index_neural_self_norm
