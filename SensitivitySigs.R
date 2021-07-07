#import dplyr for later use
library(dplyr)
library(limma)

#read GDSC data from csv file
rawData <- read.csv(file = "data2.csv")
expressionData0 <- data.frame(read.csv(file = "out_0.csv", check.names = FALSE))
expressionData1 <- data.frame(read.csv(file = "out_1.csv", check.names = FALSE))
expressionData2 <- data.frame(read.csv(file = "out_2.csv", check.names = FALSE))
expressionData3 <- data.frame(read.csv(file = "out_3.csv", check.names = FALSE))
#combine all the expression data into one frame
combinedEXP <- rbind(expressionData0, expressionData1, expressionData2, expressionData3)

#turn the raw IC50 data into a df
frame <- data.frame(rawData)


#filter frame to only include cisplatin
cis <- frame[frame$DRUG_NAME == "Cisplatin", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]

#filter frame to only include THCA, LIHC, and SKCM sample types
cisTHCA <- cis[cis$TCGA_DESC == "THCA", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]
cisLIHC <- cis[cis$TCGA_DESC == "LIHC", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]
cisSKCM <- cis[cis$TCGA_DESC == "SKCM", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]

#reassemble the dataframe
filteredCis <- rbind(cisTHCA, cisLIHC)
filteredCis <- rbind(filteredCis, cisSKCM)

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
#set up boolean explanatory variable for sensitivity
top['sensitivity'] <- TRUE
bottom['sensitivity'] <- FALSE

#recombines the top and bottom 20%
allSens <- rbind(top, bottom)
#selects only the things we need for limma
sensForJoin <- allSens[, c("COSMIC_ID", "sensitivity")]
#Do the inner join
matchedSens <- inner_join(combinedEXP, sensForJoin, by = c("GENE_SYMBOLS" = "COSMIC_ID"))
print(allSens)



#limma stuff from here on out
#make the design matrix
design <- model.matrix(~sensitivity, data = matchedSens)

#check to see that it looks right
print(head(design, 2))
#check the number of samples against the number of samples in the original set
colSums(design)
table(matchedSens[, "sensitivity"])

#remove gene symbols column
matchedSens <- matchedSens[,-1]
#remove sensitivity column
matchedSens <- matchedSens[,1:(ncol(matchedSens)-1)]

#use limma to do DE
fit <- lmFit(t(matchedSens[,-1]), design)
fit <- eBayes(fit)

#This code is to rotate the dataset later if I need

#transposed <- as.data.frame(t(matchedSens[,-1]))
#colnames(transposed) <- rownames(matchedSens)
#rownames(transposed) <- colnames(matchedSens)