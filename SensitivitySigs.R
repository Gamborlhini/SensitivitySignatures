#import dplyr, limma, samr, multtest for later use
library(dplyr)
library(limma)
library(samr)
library(multtest)

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
filteredCis <- cis[cis$TCGA_DESC == "THCA" | cis$TCGA_DESC == "LIHC" | cis$TCGA_DESC == "SKCM", c("CELL_LINE_NAME", "COSMIC_ID", "TCGA_DESC", "DRUG_NAME","LN_IC50")]


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



#limma stuff from here on out
#make the design matrix
design <- model.matrix(~sensitivity, data = matchedSens)

#rename design variable
colnames(design) <- c("sens", "sensvsres")

#remove gene symbols column
matchedSens <- matchedSens[,-1]
#remove sensitivity column
matchedSens <- matchedSens[,1:(ncol(matchedSens)-1)]

#use limma to do DE
fit <- lmFit(t(matchedSens[,-1]), design)
fit <- eBayes(fit)

#summarize results
topTable(fit, number = Inf, lfc=0.2, p.value=0.05, adjust.method="none", coef="sensvsres")

#store significantly up/downregulated genes from limma
limmaSigUp <- fit$t[fit$t[,"sensvsres"] > 1.96,]
limmaSigDown <- fit$t[fit$t[,"sensvsres"] < -1.96,]

#multtest stuff

#get the equivalent of a design variable for multtest
exp.cl <- design[,2]
#get the expression set for multtest
exp <- t(matchedSens[,-1])

#grab the stats
resP<-mt.minP(exp,exp.cl)

#store sig up/down regulated genes from multtest
multSigUp <- resP[resP$teststat > 1.96,]
multSigDown <- resP[resP$teststat < -1.96,]

#order limma up/down regulated genes
limmaSigUp <- limmaSigUp[order(-limmaSigUp[,"sensvsres"]),]
limmaSigDown <- limmaSigDown[order(limmaSigUp[,"sensvsres"]),]

#order multtest up/down regulated genes
multSigUp <- multSigUp[order(-multSigUp$teststat),]
multSigDown <- multSigDown[order(multSigUp$teststat),]

#samr stuff
exp.cl <- exp.cl + 1

samrExp <- list(x=exp, y=exp.cl, geneid=as.character(1:nrow(exp)))
samrStats <- samr(samrExp, resp.type = "Two class unpaired", assay.type = "seq")


