#Data preparation from PATRIC files

#Libraries
library(tidyverse)
#source("scripts/functions_CR.R")
library(collatRal)

####Import data####
#Load abbreviations file
meta_antibio_abbr <- read.csv("data/meta_antibio_abbr.csv", header = T, stringsAsFactors = F)

#Pick species
species <- "Escherichia-coli"
#Choose minimal number of MIC measurements to include drug
na_remove <- 200

#Load the raw data (long format) 
raw_MIC <- read.table(paste0("data/raw/", species ,".txt"), header = T, sep = "\t", dec = ",", stringsAsFactors = F)
colnames(raw_MIC) <- sapply(strsplit(colnames(raw_MIC), "\\."), function(x) x[2])
colnames(raw_MIC)[colnames(raw_MIC) == "measurement_value"] <- "MIC"

####Data cleaning####
#Remove all rows without MIC value
raw_MIC <- raw_MIC[raw_MIC$MIC != "", ]

#Create clean MIC values data frame
MIC_df <- raw_MIC %>%
#clean strings
  mutate(antibiotic = tolower(str_replace_all(antibiotic, "-", "/")),
         genome_id = str_replace_all(genome_id, "\\.", "_")) %>% 
#add antibiotic Abbreviations
  left_join(meta_antibio_abbr, by = "antibiotic") %>% 
  select(genome_id, genome_name, antibiotic, abbreviation, measurement_sign, MIC) %>% 
  filter(!is.na(MIC)) %>% 
#take MIC value from antibiotic (not beta lactamase blocker)
  mutate(MIC = sapply(strsplit(MIC, "/"), function(x) as.numeric(x[1])),
         measurement_sign = ifelse(measurement_sign %in% c(">", ">=", "=>"), ">", "<="),
         key = paste(genome_name, genome_id, abbreviation, sep = "_")) %>% 
  mutate(logMIC = log2(MIC)) %>% 
#remove duplicate rows
  distinct(MIC, key, .keep_all = TRUE)

missing_abb <- is.na(MIC_df$abbreviation)
if (any(missing_abb)) {
  warning(paste("For antibiotics:", paste(unique(MIC_df$antibiotic[!MIC_df$antibiotic %in% meta_antibio_abbr$antibiotic]), collapse = ", "), 
                "no abbreviation was found in the abbreviations file. In order to keep these in the analysis, add an abbreviation to the abbreviations file!",
                "Now, MICs without abbreviation were removed (removing", sum(missing_abb), "MIC values)"))
  MIC_df <- MIC_df[!missing_abb, ]
}


#Check if there are multiple measurements of a strain/antibiotic combination
MIC_df <- RemoveDuplicateMICs(MIC_df)

#Create wide format of MIC values
MIC_table <- MIC_df %>% 
  select(MIC, abbreviation, genome_id) %>% 
  spread(abbreviation, MIC) %>% 
  column_to_rownames(var = "genome_id")

#Remove antibiotics with to many na's
n_antibiotic <- apply(MIC_table, 2, function(x) sum(!is.na(x)))
MIC_clean <- MIC_table %>% 
  select(which(n_antibiotic > na_remove)) %>% 
  filter(rowSums(!is.na(.)) > 1)
 # filter_all(any_vars(!is.na(.)))

#Save cleaned data for further use
save(MIC_clean, file = paste0("data/clean/MIC_clean_", species,".Rdata"))
dim(MIC_clean)


####Data exploration####
#Create data frame with number of observations
tab_mic_AB <- as.data.frame(table(raw_MIC$antibiotic))
#Plot number of observations per antibiotic
p <- ggplot(data = tab_mic_AB, aes(x = Var1, y = Freq)) +
  geom_hline(yintercept = na_remove, colour = "light pink") +
  geom_bar(stat = "identity") +
  theme_bw()
p + coord_flip()+ scale_y_continuous(expand = c(0, 0), limits = c(0, length(unique(raw_MIC$genome_id))))
