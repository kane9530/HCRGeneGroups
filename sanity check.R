#sanity check


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