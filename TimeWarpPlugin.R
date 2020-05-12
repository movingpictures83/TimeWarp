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
}

run <- function() {


#t1 <- read.table("ViralChallenge_training_EXPRESSION_RMA.tsv", sep = "\t", header =FALSE, stringsAsFactors=FALSE)
#t2 <- read.table("ViralChallenge_training_CLINICAL.tsv", sep="\t", header = TRUE,  stringsAsFactors=FALSE)

t2[t2$STUDYID == "DEE4X H1N1",]="H1N1"
t2[t2$STUDYID == "DEE3 H1N1",]="H1N1"
t2[t2$STUDYID == "DEE2 H3N2",]="H3N2"
t2[t2$STUDYID == "DEE5 H3N2",]="H3N2"
t2[t2$STUDYID == "Rhinovirus Duke",]="Rhinovirus"
t2[t2$STUDYID == "Rhinovirus UVA",]="Rhinovirus"

t1 <<- as.data.frame(t(t1), stringsAsFactors=FALSE)
colnames(t1)[1] <<- "CEL"
# rownames(t1) <- substring(rownames(t1), 2, length(rownames(t1)))
# write.table(t1, file = "shortRMA.csv", sep = ",", row.names = FALSE, col.names = TRUE)
x <<- as.data.frame(merge(t1, t2, by ="CEL", stringsAsFactors=FALSE))

studyID <<- unique(x$STUDYID)
}

output <- function(outputfile) {
for(virus in studyID){
  v1 <- x[x[,"STUDYID"]==virus,]
  maxAc = 0
  for(k in c(2,3,4)){
    #P = read.csv(paste("cluster_",virus,"_","k",k,".csv", sep=""),sep = ",",header=TRUE, row.names = 1, stringsAsFactors=FALSE)
    P = read.csv(paste(prefix,virus,"_","k",k,".csv", sep=""),sep = ",",header=TRUE, row.names = 1, stringsAsFactors=FALSE)
    avgAc = 0
    for(i in k){
      lis = P[i,!is.na(P[i,])]
      grp = v1[v1$SUBJECTID %in% lis,]
      times = unique(grp$TIMEHOURS)
      grpAc = 0
      if(min(grp$SYMPTOMATIC_SC2)!= max(grp$SYMPTOMATIC_SC2))
      {
        for(t in times)
        {
          X = data.matrix(grp[grp$TIMEHOURS==t,2:22278])
          Y = as.factor(grp[grp$TIMEHOURS==t,22286])
          my_pls1 = plsda(X, Y, ncomp = 2)
          res = order(varImp(my_pls1), decreasing = TRUE)[1:50]
          
          datX = X[,c(as.numeric(unlist(res)))]
          rf.label = Y
          
          cv.folds <- createMultiFolds(rf.label, k=10, times = 10)
          fit  = trainControl(method = "repeatedcv", number = 10, repeats = 10, index = cv.folds)
          res = train(x = datX, y = rf.label, method = "rf", tuneLength = 3, ntree = 1000, trControl = fit)
          
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

