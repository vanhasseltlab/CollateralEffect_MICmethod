#Data preparation from PATRIC files

#Libraries
library(tidyverse)

####Import data####
#Load abbreviations file
meta_antibio_abbr <- read.csv("data/clean/meta_antibio_abbr.csv", header = T, stringsAsFactors = F)

#Pick species
species <- "Staphylococcus-aureus"

#Load the raw data (long format) 
raw_MIC <- read.table(paste0("data/raw/",species ,".txt"), header = T, sep = "\t", dec = ",", stringsAsFactors = F)
colnames(raw_MIC) <- c("genome_name", "genome_id", "antibiotic", "MIC", "measurement_value","laboratory_typing_method","laboratory_typing_version")

#Remove all rows without MIC value
raw_MIC <- raw_MIC[raw_MIC$MIC != "", ]

####Data cleaning####
#Pick options
na_remove <- 100 #minimal number of MIC measurements for testing

#Create clean MIC values data frame
MIC_df <- raw_MIC %>%
#clean strings
  mutate(antibiotic = tolower(str_replace_all(antibiotic, "-", "/")),
         genome_id = str_replace_all(genome_id, "\\.", "_")) %>% 
#add antibiotic Abbreviations
  left_join(meta_antibio_abbr, by = "antibiotic") %>% 
  select(genome_id, genome_name, antibiotic, abbreviation, MIC) %>% 
  filter(!is.na(MIC)) %>% 
#take MIC value from antibiotic (not beta lactamase blocker)
  mutate(MIC = sapply(strsplit(MIC, "/"), function(x) as.numeric(x[1]))) %>% 
  mutate(logMIC = log2(MIC)) %>% 
#remove duplicate rows
  distinct(MIC, genome_id, abbreviation, .keep_all = TRUE)

#Create wide format of MIC values
MIC_table <- MIC_df %>% 
  select(MIC, abbreviation, genome_id) %>% 
  spread(abbreviation, MIC) %>% 
  column_to_rownames(var = "genome_id")

#Remove antibiotics with to many na's
n_antibiotic <- apply(MIC_table, 2, function(x) sum(!is.na(x)))
MIC_clean <- MIC_table %>% 
  select(which(n_antibiotic > na_remove))

#Save cleaned data for further use
save(MIC_clean, file = paste0("data/clean/MIC_clean_", species,".Rdata"))

####Data exploration####
#Create data frame with number of observations
tab_mic_AB <- as.data.frame(table(raw_MIC$antibiotic))
#Plot number of observations per antibiotic
p <- ggplot(data = tab_mic_AB, aes(x = Var1, y = Freq)) +
  geom_hline(yintercept = na_remove, colour = "light pink") +
  geom_bar(stat = "identity") +
  theme_bw()
p + coord_flip()+ scale_y_continuous(expand = c(0, 0), limits = c(0, length(unique(raw_MIC$genome_id))))
  


unique(raw_MIC$antibiotic[grepl("/", raw_MIC$MIC)])
table(raw_MIC$MIC)

data_files <- list.files("data/raw")
for (i in 1:length(data_files)) {
  raw_MIC <- read.table(paste0("data/raw/", data_files[i]), header = T, sep = "\t", dec = ",", stringsAsFactors = F)
  colnames(raw_MIC) <- c("genome_name", "genome_id", "antibiotic", "MIC", "measurement_value","laboratory_typing_method","laboratory_typing_version")
  print(data_files[i])
  print("number of genome ids:")
  print((table(raw_MIC$antibiotic)))

}



