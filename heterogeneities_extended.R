
#Importing libraries

library(plyr)
library(tidyverse)
library(stringr)
library(qdapRegex)
library(viridis)

#Method

processCsv <- function(myDf, geneNames, mySampleNo){
  #ProcessCsv takes 3 arguments:
  # 1) myCsvFiles is the csv file exported from imaris, in the data.frame format.
  # 2) geneNames is a vector containing the gene names.
  # 3) mySampleNo is an integer indicating the sample number.
  # And returns a dataframe that has the relevant columns, with modified column names
  myDf$sampleID = mySampleNo
  csv <- select(myDf, starts_with("Intensity"), "ID", "sampleID")
  colnames(csv)[1] <- geneNames
  colnames(csv)[2] <- "cellID"
  csv <- as.data.frame(csv)
  return (csv)
}



#Listing the directories for 18ss

dir18ss <- list.dirs(path = "./dataCsv/18ss")[-1]

#Listing the full names of all csv file names within specified directories
myDirCsv<- dir18ss %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 


#Extracting the sample ID from the file names
mySampleNo <- c()
for (dir18ss in myDirCsv){
  myID <- rm_between(dir18ss, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  mySampleNo <- c(mySampleNo, myID)
}
mySampleNo <- as.integer(mySampleNo)

# Specify the gene names
geneNames = c("sox2", "tbxta")


# Obtain a list containing the csvs. Two lapply are used because we have to extract the individual elements of the list to access the csv files.
# See https://stackoverflow.com/questions/1169456/the-difference-between-bracket-and-double-bracket-for-accessing-the-el for subsetting details.

myCsvFiles <- lapply(myDirCsv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 



# Processing the individual csv files 

# Initialising an empty vector to contain all the processed Csvs
myCsvFilesProc <- vector("list", length = length(myCsvFiles))
myCsvFilesProc_merged <- vector("list", length = length(myCsvFiles))

# Looping through all elements of list

for (i in 1:length(myCsvFiles)){
  for (j in 1:length(geneNames)){
    myCsvProcessed <- processCsv(myDf = myCsvFiles[[i]][[j]], 
                                 geneNames = geneNames[j],
                                 mySampleNo = mySampleNo[i])
    myCsvFilesProc[[i]][[j]] <- myCsvProcessed
  }
  myCsvFilesProc_merged[[i]] <- merge(myCsvFilesProc[[i]][[1]], myCsvFilesProc[[i]][[2]], by = c("cellID", "sampleID"))
}



#Combine all csv into one tidied csv. Also, min-max normalisation for each column and NM index calculation.
#Index defined to be difference between sox2 and tbxta expression: negative values indicates higher tbxta to sox2 exp.

df_combined <- myCsvFilesProc_merged %>%
  reduce(rbind) %>%
  group_by(sampleID) %>%
  arrange(sampleID, cellID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nmIndex = .[[5]]-.[[6]])

df_combined

#Plotting the distribution

p <- ggplot(df_combined, aes(x=nmIndex, fill = factor(sampleID))) + geom_histogram(binwidth=0.2, color = 'black') 
p1 <- p + scale_x_continuous(name="NMindex", limits=c(-1, 1))
p2 <- p1 +scale_fill_manual(values=plasma(n=9))
p2

#time to extend this

#Listing the directories for 21ss
dir21ss <- list.dirs(path = "./dataCsv/21ss")[-1]

myDir21ssCsv<- dir21ss %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

mySampleNo <- c()
for (dir21ss in myDir21ssCsv){
  myID <- rm_between(dir21ss, "sample_", "/", extract=TRUE) %>%
    unlist() %>%
    unique()
  mySampleNo <- c(mySampleNo, myID)
}
mySampleNo <- as.integer(mySampleNo)

geneNames = c("sox2", "tbxta")

myCsvFiles <- lapply(myDir21ssCsv, function(x){ lapply(x, FUN = read.csv, header = TRUE, skip = 3)}) 

myCsvFilesProc <- vector("list", length = length(myCsvFiles))
myCsvFilesProc_merged <- vector("list", length = length(myCsvFiles))


for (i in 1:length(myCsvFiles)){
  for (j in 1:length(geneNames)){
    myCsvProcessed <- processCsv(myDf = myCsvFiles[[i]][[j]], 
                                 geneNames = geneNames[j],
                                 mySampleNo = mySampleNo[i])
    myCsvFilesProc[[i]][[j]] <- myCsvProcessed
  }
  myCsvFilesProc_merged[[i]] <- merge(myCsvFilesProc[[i]][[1]], myCsvFilesProc[[i]][[2]], by = c("cellID", "sampleID"))
}


df_combined <- myCsvFilesProc_merged %>%
  reduce(rbind) %>%
  group_by(sampleID) %>%
  arrange(sampleID, cellID) %>%
  mutate_at(.funs = list(normalise = ~((.-min(.))/max(.-min(.)))), .vars = 3:4) %>%
  ungroup()%>%
  mutate(nmIndex = .[[5]]-.[[6]])

df_combined

p <- ggplot(df_combined, aes(x=nmIndex, fill = factor(sampleID))) + geom_histogram(binwidth=0.2, color = 'black') 
p1 <- p + scale_x_continuous(name="NMindex", limits=c(-1, 1))
p2 <- p1 +scale_fill_manual(values=plasma(n=12))
p2

#the question is, can I automate this.
#he says use a lapply function
dir21ss <- list.dirs(path = "./dataCsv/21ss")[-1]

myDir21ssCsv<- dir21ss %>% 
  map(~list.files(path = ., pattern="\\.csv$", full.names = TRUE)) 

somitestages <- list.dirs(path="./dataCsv", recursive = FALSE)

what <- somitestages %>% 
  map(~list.files(path = ., pattern = "\\.csv$", full.names = TRUE)

for (somitestage in datacsvfile)