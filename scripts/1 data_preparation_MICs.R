#Data preparation from PATRIC files

#Libraries
library(tidyverse)

####Import data####
#Load the raw data (long format) 
raw_MIC <- read.table("data/raw/Klebsiella-pneumoniae.txt", header = T, sep = "\t", dec = ",", stringsAsFactors = F)
colnames(raw_MIC) <- c("genome_name", "genome_id", "antibiotic", "MIC", "measurement_value","laboratory_typing_method","laboratory_typing_version")

#Load abbreviations file
meta_antibio_abbr <- read.csv("data/clean/meta_antibio_abbr.csv", header = T)

####Data cleaning####
#Pick options
na_remove <- 100 #minimal number of MIC measurements for testing

#Create clean MIC values data frame
MIC_df <- raw_MIC %>%
  left_join(meta_antibio_abbr) %>% # add antibiotic Abbreviations
  select(genome_id, genome_name, antibiotic, abbreviation, MIC) %>% 
  filter(!is.na(MIC)) %>% 
  filter(!(grepl("/", antibiotic) | grepl("-", antibiotic))) %>% #Remove combination therapy! 
  mutate(MIC = as.numeric(MIC)) %>% 
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
save(MIC_clean, file = "data/clean/MIC_clean.Rdata")


####Data exploration####
#Create data frame with number of observations
tab_mic_AB <- as.data.frame(table(raw_MIC_rm$antibiotic))
#Plot number of observations per antibiotic
p <- ggplot(data = tab_mic_AB, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  theme_bw()
p + coord_flip()+ scale_y_continuous(expand = c(0, 0), limits = c(0, length(unique(raw_MIC_rm$genome_id))))


unique(raw_MIC$antibiotic[grepl("/", raw_MIC$MIC)])
table(raw_MIC$MIC)

