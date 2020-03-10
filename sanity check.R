#sanity check

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
