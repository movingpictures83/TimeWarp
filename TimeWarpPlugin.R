library("caret")
library("spls")

#setwd("appData/Res_Challenge/Data")
#rm(list=ls(all=TRUE))

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  print("READING INPUT FILES...");
  t1 <<- read.table(toString(parameters["training",2]), sep = "\t", header =FALSE, stringsAsFactors=FALSE)#, nrow=20000)
  t2 <<- read.table(toString(parameters["clinical",2]), sep="\t", header = TRUE,  stringsAsFactors=FALSE)
  print("DONE");
  prefix <<- toString(parameters["prefix", 2]);
  joinby <<- toString(parameters["joinby", 2])
  orderby <<- toString(parameters["orderby", 2])
  id <<- toString(parameters["id", 2])
  classcol <<- toString(parameters["classcol", 2]) 
  traincontrolmethod <<- toString(parameters["trainControl", 2]) 
  trainmethod <<- toString(parameters["train", 2]) 
  myX <<- toString(parameters["x", 2])
}

run <- function() {


#t1 <- read.table("ViralChallenge_training_EXPRESSION_RMA.tsv", sep = "\t", header =FALSE, stringsAsFactors=FALSE)
#t2 <- read.table("ViralChallenge_training_CLINICAL.tsv", sep="\t", header = TRUE,  stringsAsFactors=FALSE)

t1 <<- as.data.frame(t(t1), stringsAsFactors=FALSE)
colnames(t1)[1] <<- joinby
# rownames(t1) <- substring(rownames(t1), 2, length(rownames(t1)))
# write.table(t1, file = "shortRMA.csv", sep = ",", row.names = FALSE, col.names = TRUE)
x <<- as.data.frame(merge(t1, t2, by =joinby, stringsAsFactors=FALSE))
studyID <<- unique(as.character(unlist(x[id])))
   train_set_size <<- ncol(t1)
   class_index <<- grep(classcol, colnames(t2))


#studyID <<- unique(x$STUDYID)
}

output <- function(outputfile) {
for(virus in studyID){
  v1 <- x[x[,id]==virus,]
  maxAc = 0
  for(k in c(2,3,4)){
    P = read.csv(paste(prefix,virus,"_","k",k,".csv", sep=""),sep = ",",header=TRUE, row.names = 1, stringsAsFactors=FALSE)
    avgAc = 0
    for(i in k){
      lis = P[i,!is.na(P[i,])]
      grp = v1[as.character(unlist(v1[orderby])) %in% lis,]
      times = unique(as.character(unlist(grp[myX])))
      grpAc = 0
      if(min(as.integer(unlist(grp[classcol]))) != max(as.integer(unlist(grp[classcol]))))
      {
        for(t in times)
        {
          X = data.matrix(grp[as.double(unlist(grp[myX]))==t,2:train_set_size])
          Y = as.factor(grp[as.double(unlist(grp[myX]))==t,(train_set_size+class_index)])
          my_pls1 = plsda(X, Y, ncomp = 2)
          res = order(varImp(my_pls1), decreasing = TRUE)[1:50]
          
          datX = X[,c(as.numeric(unlist(res)))]
          rf.label = Y
          
          cv.folds <- createMultiFolds(rf.label, k=10, times = 10)
          fit  = trainControl(method = traincontrolmethod, number = 10, repeats = 10, index = cv.folds)
          res = train(x = datX, y = rf.label, method = trainmethod, tuneLength = 3, ntree = 1000, trControl = fit)
          
          # df = data.frame(max(res$results$Accuracy), t, virus)
          write.csv(lis, file = outputfile, row.names = FALSE)
          
          if(max(res$results$Accuracy)>grpAc)
          {
            grpAc = max(res$results$Accuracy)
            bestFit = fit
            bestRes = res
            bestTime = t
          }
        }
      }
      else
      {
        grpAc = 1
      }
      avgAc = grpAc/k + avgAc
    }
    if(avgAc > maxAc)
        maxAc = avgAc
  }
  print(virus)
  print(maxAc)
}
}

