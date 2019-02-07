#===================================================================================================
# This script is to do preprocessing of GDSC data and obtain a complete dataset.
# You need load the datasets from ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-5.0/. But downloading and transforming the three used datasets below to *.csv files first.
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 07-Feb-2019
#===================================================================================================
features <- data.frame(read.csv("gdsc_en_input_w5.csv", head=T))
names.fea <- strsplit(rownames(features), "")
features <- t(features)
p <- c(13321, 13747-13321, 13818-13747)
Cell.Line <- rownames(features)
features <- data.frame(Cell.Line, features)

ic50 <- data.frame(read.csv("gdsc_drug_sensitivity_fitted_data_w5.csv", head=T))
IC50IntAucAlphaBetaIQ <- ic50[,c(1,4,7,8,9,11,12,15,54,104)]
ic50 <- ic50[,c(1,4,7,11)]
drug.id <- data.frame(read.csv("gdsc_tissue_output_w5.csv", head=T))[,c(1,3)]
drug.id2 <- drug.id[!duplicated(drug.id$drug.id),]
# delete drug.id=1066 since ID1066 and ID156 both correspond drug AZD6482, and no ID1066 in the "suppl.Data1" 
drug.id2 <- drug.id2[drug.id2$drug.id!=1066,]
drug.id2$drug.name <- as.character(drug.id2$drug.name)
drug.id2$drug.name <- substr(drug.id2$drug.name, 1, nchar(drug.id2$drug.name)-6)
drug.id2$drug.name <- gsub(" ", "-", drug.id2$drug.name)

library(plyr)
# mapping the drug_id to drug names in drug sensitivity data set
ic50$drug_id <- mapvalues(ic50$drug_id, from = drug.id2[,2], to = drug.id2[,1])
colnames(ic50)[c(1,2,4)] <- c("Cell.Line", "compound", "IC50")

# transform drug sensitivity overall cell lines to a data matrix
y0 <- reshape(ic50[,c(1,2,4)], v.names="IC50", timevar="compound", idvar="Cell.Line", direction="wide")
y0$Cell.Line <- gsub("-", ".", y0$Cell.Line)

#================
# Add other information related to IC50 to do the sensitivity analysis later
#================
names_ic_other <- c("50l", "50h", "25", "75", "Beta")
indx_ic_other <- c(5, 4, 9, 10, 8)
for(i in 1:length(names_ic_other)){
  assign(paste("ic",names_ic_other[i],sep=""), IC50IntAucAlphaBetaIQ[,c(1,2,indx_ic_other[i])])
  temp <- eval(parse(text=paste("ic",names_ic_other[i],sep="")))
  temp$drug_id <- mapvalues(temp$drug_id, from = drug.id2[,2], to = drug.id2[,1])
  colnames(temp)[1:3] <- c("Cell.Line", "compound", paste("IC",names_ic_other[i],sep=""))
  temp2 <- reshape(temp[,1:3], v.names=paste("IC",names_ic_other[i],sep=""), timevar="compound", idvar="Cell.Line", direction="wide")
  temp2$Cell.Line <- gsub("-", ".", temp2$Cell.Line)
  temp2 <- temp2[temp2$Cell.Line %in% names.cell.line, c(1,1+which(substr(colnames(temp2), 4+nchar(names_ic_other[i]), nchar(colnames(temp2)))[-1] %in% names.drug))]
  temp2 <- temp2[sort.list(temp2$Cell.Line),-1]
  assign(paste("y0",names_ic_other[i],sep=""),  temp2)
}

#===============
# select nonmissing pharmacological data
#===============
y00 <- y0
m0 <- dim(y0)[2]-1
eps <- 0.05
# r1.na is better to be not less than r2.na
r1.na <- 0.3
r2.na <- 0.2
k <- 1
while(sum(is.na(y0[,2:(1+m0)]))>0){
  r1.na <- r1.na - eps/k
  r2.na <- r1.na - eps/k
  k <- k + 1
  ## select drugs with <30% (decreasing with k) missing data overall cell lines
  na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y<r1.na)<m0){
    y0 <- y0[,-c(1+which(na.y>=r1.na))]
    m0 <- sum(na.y<r1.na)
    na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  }
  
  ## select cell lines with treatment of at least 80% (increasing with k) drugs
  na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y0<r2.na)<(dim(y0)[1])){
    y0 <- y0[na.y0<r2.na,]
    na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  }
  num.na <- sum(is.na(y0[,2:(1+m0)]))
  cat("#{NA}=", num.na, "\n", "r1.na =", r1.na, ", r2.na =", r2.na, "\n")
}

#===============
# combine drug sensitivity, tissues and molecular features
#===============
yx <- merge(y0, features, by="Cell.Line")
names.cell.line <- yx$Cell.Line
names.drug <- colnames(yx)[2:(dim(y0)[2])]
names.drug <- substr(names.drug, 6, nchar(names.drug))
# numbers of gene expression features, copy number festures and muatation features
p <- c(13321, 13747-13321, 13818-13747) 
num.nonpen <- 13
yx <- data.matrix(yx[,-1])
y <- yx[,1:(dim(y0)[2]-1)]
x <- cbind(yx[,dim(y0)[2]-1+sum(p)+1:num.nonpen], yx[,dim(y0)[2]-1+1:sum(p)])

## preselect gene expression explaining 50% variations over all samples
var.x1 <- apply(log(x[,num.nonpen+1:p[1]]), 2, var)
var.sort <- sort(var.x1, decreasing=TRUE)
sum.a <- cumsum(var.sort)
half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]/2))[1]]
x1.cov <- x[,num.nonpen+1:p[1]][,var.x1>=half.a]

x <- cbind(x[,1:num.nonpen], x1.cov, x[,num.nonpen+p[1]+1:(p[2]+p[3])])
p[1] <- dim(x1.cov)[2]
m <- dim(y)[2]

GDSC <- list(y=y, x=x, m=m, p=p, num.nonpen=num.nonpen, names.cell.line=names.cell.line, names.drug=names.drug)
save(GDSC, file="GDSC_complete.RData")

