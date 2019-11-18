#Recreate calculating the CrCs rank from raw data

#Load the data
raw_MIC <- read.table("data/Klebsiella-pneumoniae.txt", header = T, sep = "\t", dec = ",", stringsAsFactors = F)
colnames(raw_MIC) <- c("genome_name", "genome_id", "antibiotic", "MIC", "measurement_value","laboratory_typing_method","laboratory_typing_version")
meta_antibio_abbr <- read.csv("default_files/meta_antibio_abbr.csv", header = T)

raw_MIC_rm <- raw_MIC[raw_MIC$MIC != "", ]
strange_values <- names(table(raw_MIC_rm$MIC)[table(raw_MIC_rm$MIC) == 1])
raw_MIC_rm[raw_MIC_rm$MIC %in% strange_values, ]

table(raw_MIC_rm$MIC)
nrow(raw_MIC_rm)
length(unique(raw_MIC_rm$genome_id))
length(unique(raw_MIC_rm$antibiotic))
tab_mic_AB <- as.data.frame(table(raw_MIC_rm$antibiotic))

p <- ggplot(data = tab_mic_AB, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  theme_bw()
p + coord_flip()+ scale_y_continuous(expand = c(0, 0), limits = c(0, length(unique(raw_MIC_rm$genome_id))))


#Libaries
library(tidyverse)

#Pick options
na_remove <- 0.8 #minimal proportion of non missing measurements

#Data cleaning
#get clean MIC values data frame
MIC_df <- raw_MIC %>% #create MIC_list.txt
  mutate(MIC = sapply(strsplit(MIC, "/"), function(x) mean(as.numeric(x))),
         antibiotic = tolower(str_replace_all(antibiotic, "-", "/")),
         genome_id = str_replace_all(genome_id, "\\.", "_")) %>% 
  left_join(meta_antibio_abbr) %>% # add antibiotic Abbreviations
  select(genome_id, genome_name, antibiotic, abbreviation, MIC) %>% 
  filter(!is.na(MIC)) %>% 
  distinct(MIC, genome_id, abbreviation, .keep_all= TRUE)
rm(raw_MIC, meta_antibio_abbr)

#get wide format of MIC values
MIC_table <- MIC_df %>% 
  select(MIC, abbreviation, genome_id) %>% 
  spread(abbreviation, MIC) %>% 
  column_to_rownames(var = "genome_id")

#remove columns with to many na's
na_cols <- apply(MIC_table, 2, function(x) mean(is.na(x)))
meas_cols <- apply(MIC_table, 2, function(x) sum(!is.na(x)))
na_rows <- apply(MIC_table, 1, function(x) mean(is.na(x)))
hist(na_cols, breaks = 30)
hist(na_rows, breaks = 30)
sum(meas_cols < 100)
MIC_clean <- MIC_table %>% 
  select(which(na_cols < na_remove))


na_rows_updated <- apply(MIC_clean, 1, function(x) mean(is.na(x)))
hist(na_rows_updated, breaks = 30)

#Dichotomization on exp(mean(log(X)))
mean(MIC_clean$AMK, na.rm = T)
exp(mean(log(MIC_clean$AMK), na.rm = T))
sum(MIC_clean$AMK <= 8, na.rm = T)
sum(MIC_clean$AMK > 8, na.rm = T)
barplot(table(MIC_clean$AMK))
# for (i in 1:20){
#   barplot(table(MIC_clean[, i]), main = round(exp(mean(log(MIC_clean[, i]), na.rm = TRUE)), 2))
# }
save(MIC_clean, file = "test_new_measure_data/MIC_clean.Rdata")