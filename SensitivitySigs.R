library(dplyr)

rawData <- read.csv(file = "data2.csv")

frame <- data.frame(rawData)

print(names(frame))

cis <- frame[frame$DRUG_NAME == "Cisplatin", c("CELL_LINE_NAME","DRUG_NAME","LN_IC50")]

print(nrow(cis))

getTop <- function(x){
  return(top_frac(x, 0.2))
}

getBottom <- function(x){
  return(top_frac(x, -0.2))
}


top <- getTop(cis)
bottom <- getBottom(cis)

print(top)
print(bottom)
