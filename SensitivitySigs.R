#import dplyr for later use
library(dplyr)

#read GDSC data from csv file
rawData <- read.csv(file = "data2.csv")
expressionData0 <- data.frame(read.csv(file = "out_0.csv", check.names = FALSE))
expressionData1 <- data.frame(read.csv(file = "out_1.csv", check.names = FALSE))
expressionData2 <- data.frame(read.csv(file = "out_2.csv", check.names = FALSE))
expressionData3 <- data.frame(read.csv(file = "out_3.csv", check.names = FALSE))
combinedEXP <- rbind(expressionData0, expressionData1, expressionData2, expressionData3)

frame <- data.frame(rawData)
print(names(combinedEXP))

#filter frame to only include cisplatin
cis <- frame[frame$DRUG_NAME == "Cisplatin", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]

#filter frame to only include THCA, LIHC, and SKCM sample types
cisTHCA <- cis[cis$TCGA_DESC == "THCA", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]
cisLIHC <- cis[cis$TCGA_DESC == "LIHC", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]
cisSKCM <- cis[cis$TCGA_DESC == "SKCM", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]

#reassemble the dataframe
filteredCis <- rbind(cisTHCA, cisLIHC)
filteredCis <- rbind(filteredCis, cisSKCM)

#display the number of cell lines in the filtered frame
print(nrow(filteredCis))

#takes a dataframe and returns the top 20% based on the last column values (LN_IC50)
getTop <- function(x){
  return(top_frac(x, 0.2, LN_IC50))
}

#takes a dataframe and returns the bottom 20% based on the last column values (LN_IC50) 
getBottom <- function(x){
  return(top_frac(x, -0.2, LN_IC50))
}

#Creates dataframes for the top and bottom responders
top <- getTop(filteredCis)
bottom <- getBottom(filteredCis)
#set up boolean for sensitivity
top['sensitivity'] <- TRUE
bottom['sensitivity'] <- FALSE

allSens <- rbind(top, bottom)
#Do the inner join
matchedSens <- inner_join(allSens, combinedEXP, by = c("COSMIC_ID" = "GENE_SYMBOLS"))
print(allSens)

print(matchedSens)
