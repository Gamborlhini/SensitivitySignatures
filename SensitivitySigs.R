#import dplyr, limma, samr, multtest for later use
library(dplyr, mask.ok = T)
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

#use limma to get a linear model
fit <- lmFit(t(matchedSens), design)
#compute eBayes statistics
fit <- eBayes(fit)

#summarize results
topTable(fit, number = Inf, lfc=0.2, p.value=0.05, adjust.method="none", coef="sensvsres")

#store significantly up/downregulated genes from limma
limmaSigUp <- fit$t[fit$p.value[,"sensvsres"] < 0.005 & fit$t[,"sensvsres"] > 1.96,]
limmaSigDown <- fit$t[fit$p.value[,"sensvsres"] < 0.005 & fit$t[,"sensvsres"] < -1.96,]

#multtest stuff

#format the equivalent of a design variable for multtest
exp.cl <- design[,2]

#format the expression set for multtest
exp <- t(matchedSens)

#compute the statistics with multtest
resP<-mt.minP(exp,exp.cl)

#store sig up/down regulated genes from multtest
multSigUp <- resP[resP$rawp < 0.005 & resP$teststat > 1.96,]
multSigDown <- resP[resP$rawp < 0.005 & resP$teststat < -1.96,]

#order limma up/down regulated genes
limmaSigUp <- limmaSigUp[order(-limmaSigUp[,"sensvsres"]),]
limmaSigDown <- limmaSigDown[order(limmaSigDown[,"sensvsres"]),]

#order multtest up/down regulated genes
multSigUp <- multSigUp[order(-multSigUp$teststat),]
multSigDown <- multSigDown[order(multSigDown$teststat),]

#samr
#make into 1,2 instead of 0,1
exp.cl <- exp.cl+1
#get rid of the names in the vector, they sometimes cause problems
exp.cl <- unname(exp.cl)
#calculate the statistics
samrStats <- SAM(exp, exp.cl, resp.type = "Two class unpaired", genenames = rownames(exp), nperms = 1000, fdr.output = 0.05)
#extract p values for each gene
samrPValues <- samr.pvalues.from.perms(samrStats$samr.obj$tt, samrStats$samr.obj$ttstar)
#take all the genes where p < .005
samrSigGenes <- samrPValues[which(samrPValues < 0.005)]

#take all genes which were upregulated
samrFullUp <- names(samrStats$samr.obj$tt[samrStats$samr.obj$tt > 0])
#intersect with those which were significant
samrSigUp <- samrSigGenes[names(samrSigGenes) %in% samrFullUp]
#take all genes which were downregulated
samrFullDown <- names(samrStats$samr.obj$tt[samrStats$samr.obj$tt < 0])
#intersect with those which were significant
samrSigDown <- samrSigGenes[names(samrSigGenes) %in% samrFullDown]

#only the names
samrUpNames <- names(samrSigUp)
samrDownNames <- names(samrSigDown)

#only the names
multUpNames <- rownames(multSigUp)
multDownNames <- rownames(multSigDown)

#only the names
limmaUpNames <- rownames(limmaSigUp)
limmaDownNames <- rownames(limmaSigDown)

#take the intersection of all of the upregulated names
seedGenesUp <- samrUpNames[samrUpNames %in% multUpNames]
seedGenesUp <- seedGenesUp[seedGenesUp %in% limmaUpNames]

#take the intersection of all the downregulated names
seedGenesDown <- samrDownNames[samrDownNames %in% multDownNames]
seedGenesDown <- seedGenesDown[seedGenesDown %in% limmaDownNames]
