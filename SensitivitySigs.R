library(dplyr)

rawData <- read.csv(file = "data2.csv")

frame <- data.frame(rawData)

print(names(frame))

cis <- frame[frame$DRUG_NAME == "Cisplatin", c("CELL_LINE_NAME", "TCGA_DESC", "DRUG_NAME","LN_IC50")]
cis <- cis[cis$TCGA_DESC != "UNCLASSIFIED", c("CELL_LINE_NAME", "TCGA_DESC", "DRUG_NAME","LN_IC50")]

cisTHCA <- cis[cis$TCGA_DESC == "THCA", c("CELL_LINE_NAME", "TCGA_DESC", "DRUG_NAME","LN_IC50")]
cisLIHC <- cis[cis$TCGA_DESC == "LIHC", c("CELL_LINE_NAME", "TCGA_DESC", "DRUG_NAME","LN_IC50")]
cisSKCM <- cis[cis$TCGA_DESC == "SKCM", c("CELL_LINE_NAME", "TCGA_DESC", "DRUG_NAME","LN_IC50")]

filteredCis <- rbind(cisTHCA, cisLIHC)
filteredCis <- rbind(filteredCis, cisSKCM)

print(nrow(filteredCis))

getTop <- function(x){
  return(top_frac(x, 0.2))
}

getBottom <- function(x){
  return(top_frac(x, -0.2))
}


top <- getTop(filteredCis)
bottom <- getBottom(filteredCis)

print(top)
print(bottom)
