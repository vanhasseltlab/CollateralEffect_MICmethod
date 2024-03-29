#!/bash/bin

################################################################################
# Data retrieval from PATRIC for replication of "Identification of collateral  #
# effects of antibiotic resistance in population surveillance data"            #
# Authors: Yob Haakman & Laura B. Zwep                                         #
################################################################################


##################----IMPORT VARIABLES---############################################
# Rename the numeric variables to their original names, for better script readabiliy.
filename=${1}

if [[ $(dpkg-query -W patric-cli | wc -l) -eq 0 ]];then   #Checks for PCLI
  curl -O -L https://github.com/PATRIC3/PATRIC-distribution/releases/download/1.018/patric-cli-1.018.deb
  sudo dpkg -i patric-cli-1.018.deb
  sudo sudo apt-get -y -f install
  echo "\nThe Patric Command Line Interface is installed.\n"
fi


##################----DOWNLOAD MICDATA---############################################
while read org
do
  # if org has no usefull data available (wc -l -eq 0) then move on to the next org
  if [[ $(p3-drug-amr-data --equal laboratory_typing_method,MIC --equal genome_name,${org} | wc -l) -lt 500 ]];then
    echo "No drug data for ${org}"
    continue
  fi
  # Get drug AMR data, where the laboratory typing method is Minimum Inhibitory Concentration (MIC) and genome_id == org. Save to file (sorted on column 1, genome_name)
  echo "${org}"
## Fields available in p3-drug-amr-data, specify after --attr
# $0   	display all columns
# $1 	  genome_id	
# $2 	  genome_name	
# $3 	  taxon_id
# $4 	  antibiotic
# $5 	  resistant_phenotype
# $6 	  measurement
# $7 	  measurement_sign
# $8   	measurement_value	
# $9 	  measurement_unit
# $10 	laboratory_typing_method
# $11 	laboratory_typing_method_version
  p3-drug-amr-data --equal laboratory_typing_method,MIC --equal genome_name,${org} --ne laboratory_typing_platform,population\ analysis \
  --attr genome_name,genome_id,antibiotic,measurement_sign,measurement_value,measurement_unit,resistant_phenotype,laboratory_typing_platform > data/${org}.txt 
done < ${filename}


