
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
  #my_csv_files_proc_merged_18[[i]] <- bind_rows(list(my_csv_files_proc_18[i][1],my_csv_files_proc_18[i][2],my_csv_files_proc_18[i][3]))
}


df_combined_18 <- my_csv_files_proc_merged_18 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#View(df_combined_18)

p_18 <- ggplot(df_combined_18, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth=0.2, color = 'black')+
  scale_x_continuous( limits=c(-1, 1))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=length(my_csv_files_18)))+
  theme_minimal() + 
  labs(x="NM index", y="Number of Cells")+
  ggtitle("18 somite stage") + 
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))
p_18

line_18 <- ggplot(df_combined_18, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+ 
  ggtitle("18 somite stage")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")
line_18


# 21ss
processcsv21 <- function(my_df, gene_names, my_sample_no_21){
  my_df$sample_ID = my_sample_no_21
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

dir21ss <- list.dirs(path = "./dataCsv/21ss")[-1]

my_dir21ss_csv<- dir21ss %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_no_21 <- c()
for (dir21ss in my_dir21ss_csv){
  my_ID_21 <- rm_between(dir21ss, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_no_21 <- c(my_sample_no_21, my_ID_21)
}

my_sample_no_21 <- as.integer(my_sample_no_21)

gene_names = c("sox2", "tbxta")

my_csv_files_21 <- lapply(my_dir21ss_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_21 <- vector("list", length = length(my_csv_files_21))
my_csv_files_proc_merged_21 <- vector("list", length = length(my_csv_files_21))

for (i in 1:length(my_csv_files_21)){
  for (j in 1:length(gene_names)){
    my_csv_processed_21 <- processcsv21(my_df = my_csv_files_21[[i]][[j+1]], 
                                        gene_names = gene_names[j],
                                        my_sample_no_21 = my_sample_no_21[i])
    my_csv_files_proc_21[[i]][[j]] <- my_csv_processed_21
  }
  my_csv_files_proc_merged_21[[i]] <- merge(my_csv_files_proc_21[[i]][[1]], my_csv_files_proc_21[[i]][[2]], by = c("cell_ID", "sample_ID"))
  #my_csv_files_proc_merged_21[[i]] <- bind_rows(list(my_csv_files_proc_21[i][1],my_csv_files_proc_21[i][2],my_csv_files_proc_21[i][3]))
}


df_combined_21 <- my_csv_files_proc_merged_21 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#View(df_combined_21)

p_21 <- ggplot(df_combined_21, aes(x=nm_index, fill = factor(sample_ID))) + geom_histogram(binwidth=0.2, color = 'black') 
p1_21 <- p_21 + scale_x_continuous(name="NMindex", limits=c(-1, 1)) + theme_minimal()
p2_21 <- p1_21 +scale_fill_manual(values=plasma(n=length(my_csv_files_21)))+theme(text=element_text(size = 15, family = "sans"))
p3_21 <- p2_21 + ggtitle("21 somite stage") + theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
p3_21

line_21 <- ggplot(df_combined_21, aes(x=nm_index))+
  geom_density(colour = "orange")+ theme_minimal()+ ggtitle("21 somite stage")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size = 15, family = "sans"))
line_21

# 24ss 

processcsv24 <- function(my_df, gene_names, my_sample_no_24){
  my_df$sample_ID = my_sample_no_24
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

dir24ss <- list.dirs(path = "./dataCsv/24ss")[-1]

my_dir24ss_csv<- dir24ss %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 


my_sample_no_24 <- c()
for (dir24ss in my_dir24ss_csv){
  my_ID_24 <- rm_between(dir24ss, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_no_24 <- c(my_sample_no_24, my_ID_24)
}

my_sample_no_24 <- as.integer(my_sample_no_24)

gene_names = c("sox2", "tbxta")

my_csv_files_24 <- lapply(my_dir24ss_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_24 <- vector("list", length = length(my_csv_files_24))
my_csv_files_proc_merged_24 <- vector("list", length = length(my_csv_files_24))

for (i in 1:length(my_csv_files_24)){
  for (j in 1:length(gene_names)){
    my_csv_processed_24 <- processcsv24(my_df = my_csv_files_24[[i]][[j+1]], 
                                        gene_names = gene_names[j],
                                        my_sample_no_24 = my_sample_no_24[i])
    my_csv_files_proc_24[[i]][[j]] <- my_csv_processed_24
  }
  my_csv_files_proc_merged_24[[i]] <- merge(my_csv_files_proc_24[[i]][[1]], my_csv_files_proc_24[[i]][[2]], by = c("cell_ID", "sample_ID"))
  #my_csv_files_proc_merged_24[[i]] <- bind_rows(list(my_csv_files_proc_24[i][1],my_csv_files_proc_24[i][2],my_csv_files_proc_24[i][3]))
}


df_combined_24 <- my_csv_files_proc_merged_24 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#View(df_combined_24)

p_24 <- ggplot(df_combined_24, aes(x=nm_index, fill = factor(sample_ID))) + geom_histogram(binwidth=0.2, color = 'black') 
p1_24 <- p_24 + scale_x_continuous(name="NMindex", limits=c(-1, 1)) + theme_minimal()
p2_24 <- p1_24 +scale_fill_manual(values=plasma(n=length(my_csv_files_24)))+theme(text=element_text(size = 15, family = "sans"))
p3_24 <- p2_24 + ggtitle("24 somite stage") + theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
p3_24

line_24 <- ggplot(df_combined_24, aes(x=nm_index))+
  geom_density(colour = "gold")+ theme_minimal()+ ggtitle("24 somite stage")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size = 15, family = "sans"))
line_24
# 26-28ss
processcsv2628ss <- function(my_df, gene_names, my_sample_no_2628ss){
  my_df$sample_ID = my_sample_no_2628ss
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

dir2628ss <- list.dirs(path = "./dataCsv/26-28ss")[-1]

my_dir2628ss_csv<- dir2628ss %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_no_2628ss <- c()
for (dir2628ss in my_dir2628ss_csv){
  my_ID_2628ss <- rm_between(dir2628ss, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_no_2628ss <- c(my_sample_no_2628ss, my_ID_2628ss)
}

my_sample_no_2628ss <- as.integer(my_sample_no_2628ss)

gene_names = c("sox2", "tbxta")

my_csv_files_2628ss <- lapply(my_dir2628ss_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_2628ss <- vector("list", length = length(my_csv_files_2628ss))
my_csv_files_proc_merged_2628ss <- vector("list", length = length(my_csv_files_2628ss))

for (i in 1:length(my_csv_files_2628ss)){
  for (j in 1:length(gene_names)){
    my_csv_processed_2628ss <- processcsv2628ss(my_df = my_csv_files_2628ss[[i]][[j+1]], 
                                                gene_names = gene_names[j],
                                                my_sample_no_2628ss = my_sample_no_2628ss[i])
    my_csv_files_proc_2628ss[[i]][[j]] <- my_csv_processed_2628ss
  }
  my_csv_files_proc_merged_2628ss[[i]] <- merge(my_csv_files_proc_2628ss[[i]][[1]], my_csv_files_proc_2628ss[[i]][[2]], by = c("cell_ID", "sample_ID"))
  #my_csv_files_proc_merged_2628ss[[i]] <- bind_rows(list(my_csv_files_proc_2628ss[i][1],my_csv_files_proc_2628ss[i][2],my_csv_files_proc_2628ss[i][3]))
}


df_combined_2628ss <- my_csv_files_proc_merged_2628ss %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#View(df_combined_2628ss)

p_2628 <- ggplot(df_combined_2628ss, aes(x=nm_index, fill = factor(sample_ID))) + geom_histogram(binwidth=0.2, color = 'black') 
p1_2628 <- p_2628 + scale_x_continuous(name="NMindex", limits=c(-1, 1)) + theme_minimal()
p2_2628 <- p1_2628 +scale_fill_manual(values=plasma(n=length(my_csv_files_2628ss)))+theme(text=element_text(size = 15, family = "sans"))
p3_2628 <- p2_2628 + ggtitle("26-28 somite stage") + theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
p3_2628

line_2628ss <- ggplot(df_combined_2628ss, aes(x=nm_index))+
  geom_density(colour = "green")+ theme_minimal()+ ggtitle("26-28 somite stage")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size = 15, family = "sans"))

line_2628ss
# 30ss

processcsv30ss <- function(my_df, gene_names, my_sample_no_30ss){
  my_df$sample_ID = my_sample_no_30ss
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

dir30ss <- list.dirs(path = "./dataCsv/30ss")[-1]

my_dir30ss_csv<- dir30ss %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

my_sample_no_30ss <- c()
for (dir30ss in my_dir30ss_csv){
  my_ID_30ss <- rm_between(dir30ss, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  my_sample_no_30ss <- c(my_sample_no_30ss, my_ID_30ss)
}

my_sample_no_30ss <- as.integer(my_sample_no_30ss)

gene_names = c("sox2", "tbxta")

my_csv_files_30ss <- lapply(my_dir30ss_csv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

my_csv_files_proc_30ss <- vector("list", length = length(my_csv_files_30ss))
my_csv_files_proc_merged_30ss <- vector("list", length = length(my_csv_files_30ss))

for (i in 1:length(my_csv_files_30ss)){
  for (j in 1:length(gene_names)){
    my_csv_processed_30ss <- processcsv30ss(my_df = my_csv_files_30ss[[i]][[j+1]], 
                                            gene_names = gene_names[j],
                                            my_sample_no_30ss = my_sample_no_30ss[i])
    my_csv_files_proc_30ss[[i]][[j]] <- my_csv_processed_30ss
  }
  my_csv_files_proc_merged_30ss[[i]] <- merge(my_csv_files_proc_30ss[[i]][[1]], my_csv_files_proc_30ss[[i]][[2]], by = c("cell_ID", "sample_ID"))
  #my_csv_files_proc_merged_30ss[[i]] <- bind_rows(list(my_csv_files_proc_30ss[i][1],my_csv_files_proc_30ss[i][2],my_csv_files_proc_30ss[i][3]))
}


df_combined_30ss <- my_csv_files_proc_merged_30ss %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#View(df_combined_30ss)

p_30 <- ggplot(df_combined_30ss, aes(x=nm_index, fill = factor(sample_ID))) + geom_histogram(binwidth=0.2, color = 'black') 
p1_30 <- p_30 + scale_x_continuous(name="NMindex", limits=c(-1, 1)) + theme_minimal()
p2_30 <- p1_30 +scale_fill_manual(values=plasma(n=length(my_csv_files_2628ss)))+theme(text=element_text(size = 15, family = "sans"))
p3_30 <- p2_30 + ggtitle("30 somite stage") + theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
p3_30

line_30ss <- ggplot(df_combined_30ss, aes(x=nm_index, y=))+
  geom_density(colour = "blue")+ theme_minimal()+ ggtitle("30 somite stage")+theme(plot.title = element_text(hjust = 0.5))+theme(text=element_text(size = 15, family = "sans"))
line_30ss

# combine dataframes

df_final_18 <- df_combined_18 %>% 
  add_column(somite_stage = somitestages[1]) %>% 
  add_column(ssnumber = 1)

df_final_18$ssnumber <- factor(df_final_18$ssnumber)

df_final_21 <- df_combined_21 %>% 
  add_column(somite_stage = somitestages[2])%>% 
  add_column(ssnumber = 2)

df_final_21$ssnumber <- factor(df_final_21$ssnumber)

df_final_24 <- df_combined_24 %>% 
  add_column(somite_stage = somitestages[3])%>% 
  add_column(ssnumber = 3)

df_final_24$ssnumber <- factor(df_final_24$ssnumber)

df_final_2628 <- df_combined_2628ss %>% 
  add_column(somite_stage = somitestages[4])%>% 
  add_column(ssnumber = 4)

df_final_2628$ssnumber <- factor(df_final_2628$ssnumber)

df_final_30 <- df_combined_30ss %>% 
  add_column(somite_stage = somitestages[5])%>% 
  add_column(ssnumber = 5)

df_final_30$ssnumber <- factor(df_final_30$ssnumber)

df_combined_final <- df_final_18 %>% 
  rbind(df_final_21) %>% 
  rbind(df_final_24) %>% 
  rbind(df_final_2628) %>% 
  rbind(df_final_30)

# combine/compare graphs 

ggarrange(p_18, p3_21 + rremove("ylab"), p3_24 + rremove("ylab"), p3_2628 + rremove("ylab"), p3_30 + rremove("ylab"),
          ncol = 5, nrow = 1)

ggarrange(line_18, line_21  + rremove("ylab"), line_24 + rremove("ylab"), line_2628ss + rremove("ylab"), line_30ss + rremove("ylab"), 
          ncol = 5, nrow = 1)

compare_lines <- ggplot(df_combined_final, aes(x=nm_index, y=, colour=ssnumber))+
  geom_density(size=2) +
  #]theme(panel.background = element_rect(fill = "black")) + 
  labs(x="NM index", y="Density", title = "Comparison across all somite stages") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_minimal()+
  theme(text=element_text(size = 15, family = "sans"))+
  scale_color_viridis_d(labels=c("18", "21", "24", "26-28", "30"), option = "D")+
  guides(colour=guide_legend(title="Somite stage"))

compare_lines

#and then something with plyrs join(), actually maybe rowbind
#but first, for each one, use addcolumn 
files_somitestages <- list.dirs(path="./dataCsv", recursive = T)
somite_s18 <- str_detect(files_somitestages, pattern = "18ss")
somitestages <- basename(files_somitestages)# %>% 
  #str_detect(pattern = "ss")
directory <- as.character()


for ( i in 1:length(somitestages)){
  directory <- list.dirs(path = "./dataCsv")
  mydirectory <- directory[i+1]%>% 
    map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE))
}

##### Trying to plot just heterogeneity######
nmp_sox2 <- ggplot(data = df_combined_final, mapping = aes(x=sox2))+
  geom_density()

nmp_sox2

nmp_sox2_norm <- ggplot(data = df_combined_final, mapping = aes(x=sox2_normalise))+
  geom_density()

nmp_sox2_norm

nmp_tbxta <- ggplot(data = df_combined_final, mapping = aes(x=tbxta))+
  geom_density()

nmp_tbxta

nmp_tbxta_norm <- ggplot(data = df_combined_final, mapping = aes(x=tbxta_normalise))+
  geom_density()

nmp_tbxta_norm

nmp_sox2_tbxta <- ggplot(data = df_combined_final, mapping = aes())+
  geom_density(aes(x=sox2), colour="blue")+
  geom_density(aes(x=tbxta), colour="red")

nmp_sox2_tbxta

nmp_sox2_tbxta_norm <- ggplot(data = df_combined_final)+
  geom_density(aes(x=sox2_normalise, colour = "Sox2"), size = 3)+
  geom_density(aes(x=tbxta_normalise, colour = "Tbxta"), size = 3)+
  xlab("Normalised intensity values")+
  scale_colour_manual("", breaks = c("Sox2", "Tbxta"),
                      values = c("red", "blue"))+
  theme_minimal()
nmp_sox2_tbxta_norm

##### Trying to justify the heterogeneity of NMPs vs noto and neural tube ######
##notochord processing function
processcsvnoto <- function(my_df, gene_names, my_sample_no_noto){
  my_df$sample_ID = my_sample_no_noto
  csv <- select(my_df, starts_with("Intensity"), "ID", "sample_ID")
  colnames(csv)[1] <- gene_names
  colnames(csv)[2] <- "cell_ID"
  csv <- as.data.frame(csv)
  return (csv)
}

##import the notochord data

dirnotoss <- list.dirs(path = "./dataCsv/18ss")[-1]

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
  #my_csv_files_proc_merged_18[[i]] <- bind_rows(list(my_csv_files_proc_18[i][1],my_csv_files_proc_18[i][2],my_csv_files_proc_18[i][3]))
}


df_combined_18 <- my_csv_files_proc_merged_18 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#View(df_combined_18)

p_18 <- ggplot(df_combined_18, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth=0.2, color = 'black')+
  scale_x_continuous( limits=c(-1, 1))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=length(my_csv_files_18)))+
  theme_minimal() + 
  labs(x="NM index", y="Number of Cells")+
  ggtitle("18 somite stage") + 
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))
p_18

line_18 <- ggplot(df_combined_18, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+ 
  ggtitle("18 somite stage")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")
line_18
##import the neural tube data
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
  #my_csv_files_proc_merged_18[[i]] <- bind_rows(list(my_csv_files_proc_18[i][1],my_csv_files_proc_18[i][2],my_csv_files_proc_18[i][3]))
}


df_combined_18 <- my_csv_files_proc_merged_18 %>%
  reduce(rbind) %>%
  group_by(sample_ID) %>%
  arrange(sample_ID, cell_ID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nm_index = .[[5]]-.[[6]])

#View(df_combined_18)

p_18 <- ggplot(df_combined_18, aes(x=nm_index, fill = factor(sample_ID)))+
  geom_histogram(binwidth=0.2, color = 'black')+
  scale_x_continuous( limits=c(-1, 1))+
  theme(text=element_text(size = 30, family = "sans"))+
  scale_fill_manual(values=viridis(n=length(my_csv_files_18)))+
  theme_minimal() + 
  labs(x="NM index", y="Number of Cells")+
  ggtitle("18 somite stage") + 
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=guide_legend(title="Sample no."))
p_18

line_18 <- ggplot(df_combined_18, aes(x=nm_index))+
  geom_density(colour = "darkslateblue") +
  theme_classic2()+ 
  ggtitle("18 somite stage")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size = 13, family = "sans"))+
  labs(x="NM index", y="Density")
line_18