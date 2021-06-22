#import dplyr for later use
library(dplyr)


#read GDSC data from csv file
rawData <- read.csv(file = "data2.csv")
expressionData0 <- read.csv(file = "CellLineBasalExpression_0.csv")
expressionData1 <- read.csv(file = "CellLineBasalExpression_1.csv")
expressionData2 <- read.csv(file = "CellLineBasalExpression_2.csv")
expressionData3 <- read.csv(file = "CellLineBasalExpression_3.csv")
combinedEXP <- rbind(expressionData0, expressionData1, expressionData2, expressionData3)
#reformat so cell lines are rows and genes are columns
#expressionData.t <- t(expressionData)
#convert into dataframe
frame <- data.frame(rawData)
#display all headers
print(names(frame))
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
  return(top_frac(x, 0.2))
}

#takes a dataframe and returns the bottom 20% based on the last column values (LN_IC50) 
getBottom <- function(x){
  return(top_frac(x, -0.2))
}

#Creates dataframes for the top and bottom responders
top <- getTop(filteredCis)
bottom <- getBottom(filteredCis)
#set up boolean for sensitivity
top['sensitivity'] <- TRUE
bottom['sensitivity'] <- FALSE

allSens <- rbind(top, bottom)
#Displays top and bottom respondents
print(allSens)
print(expressionData0 %>% slice(1:2))
