library(caret)
library(glmnet)

###
###Get data for first fold
fold1.df<-read.table("fold1_signiSNPs_fullsamples.raw", sep=" ", header=TRUE, as.is=TRUE)
dim(fold1.df)
#get sample IDs for fold1
fold1_trainIDs<-read.table("../../input/fold1_keep.txt", sep="\t", header=FALSE, as.is=TRUE)
#there are still some missing genotypes that need to be removed so that glm works
fold1_noNAs.df<-fold1.df[!unlist(vapply(fold1.df, anyNA, logical(1)))]
dim(fold1_noNAs.df)

#get train and test using the original fold sample IDs
train1.df<-fold1_noNAs.df[fold1_noNAs.df$IID %in% fold1_trainIDs$V2,]
#to get the oposit of in, I first create a function
'%ni%' <- Negate('%in%')
test1.df<-fold1_noNAs.df[fold1_noNAs.df$IID %ni% fold1_trainIDs$V2,]

###
###GLM with LASSO penalization for feature selection
y1.train <-factor(train1.df$PHENOTYPE)
x1.train<-data.matrix(train1.df[-1:-6])
fit1 <- glmnet(x1.train,y1.train,family="binomial",standardize=TRUE)
coefs1=coef(fit1,s=min(fit1$lambda))
lasso.coeff.best1<-cbind(coefs1@Dimnames[[1]][which(coefs1 !=0)],coefs1[which(coefs1 != 0)])
head(lasso.coeff.best1)
#     [,1]           [,2]                  
#[1,] "(Intercept)"  "2.33029016149315"    
#[2,] "rs4648612_T"  "0.0889067319101327"  
#[3,] "rs4366263_G"  "1.86266120775391e-15"
#[4,] "rs28778974_T" "2.27153805823647e-16"
#[5,] "rs12738895_T" "0.0620313426477561"  
#[6,] "rs72647619_C" "-0.0113019269147701" 

x1_lasso.train<-train1.df[,colnames(train1.df) %in% lasso.coeff.best1[-1,1]]

#We now need a dependent variable of characters otherwise it will throw an error Error: At least one of the class levels is not a valid R variable name; 
#This will cause errors when class probabilities are generated because the variables names will be converted to  X1, X2 . Please use factor levels that can be used as valid R variable names  (see ?#make.names for help).
y1.train.cat <- sapply(y1.train, switch,
"1" = "No",
"2" = "Yes")
#Additionally we have to have everything in one dataframe to feed it so I merge x and y again
#> dat1.train <-cbind(x1_lasso.train,y1.train.cat)
dat1.train<-x1_lasso.train
dat1.train["PHENOTYPE"]<-y1.train.cat


set.seed(3233)
fitControl<-trainControl(method="cv",number=10, classProbs=TRUE, summaryFunction=twoClassSummary, savePredictions=TRUE)

svmFit1 <- train(PHENOTYPE ~ ., data =dat1.train, method="svmLinear2",tuneLength=10, trControl=fitControl,preProc = c("center","scale"),metric= "ROC")

#There was an observation with a missing Phenotype in the test set, so I remove it
test1.df<-test1.df[-74,]
y1.test<-factor(test1.df$PHENOTYPE)
x1_lasso.test<-test1.df[,colnames(test1.df) %in% lasso.coeff.best1[-1,1]]
dat1.test<-x1_lasso.test
y1.test.cat <- sapply(y1.test, switch,
"1" = "No",
"2" = "Yes")
dat1.test["PHENOTYPE"]<-y1.test.cat

###
###Support Vector Machine
#We will use linear and non-linear SVM
svm.pred1<-predict(svmFit1,newdata=dat1.test)
confusionMatrix(svm.pred1,as.factor(dat1.test$PHENOTYPE),positive = "YES")
       
#Try gridsearch on Cost
grid <- expand.grid(C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5))
set.seed(3233)
svm_Linear_Grid <- train(PHENOTYPE ~., data = dat1.train, method = "svmLinear",#change to svmLinear because svmLinear2 throws an error
                    trControl=fitControl,
                    preProcess = c("center", "scale"),
                    tuneGrid = grid,
                    tuneLength = 10)
 
svm_Linear_Grid
test_pred_grid <- predict(svm_Linear_Grid, newdata = dat1.test)
confusionMatrix(test_pred_grid,as.factor(dat1.test$PHENOTYPE),positive = "YES")

#Try non-lineal RBF kernel
svmFit1_RBF <- train(PHENOTYPE ~ ., data =dat1.train, method="svmRadial",tuneLength=10, trControl=fitControl,preProc = c("center","scale"),metric= "ROC")
svmFit1_RBF
test_pred_RBF <- predict(svmFit1_RBF, newdata = dat1.test)
confusionMatrix(test_pred_RBF,as.factor(dat1.test$PHENOTYPE),positive = "YES")


###
###Random Forest
#We will use the popular Random Forest algorithm as the subject of our algorithm tuning.
#Random Forest is not necessarily the best algorithm for this dataset, but it is a very popular algorithm
#In this case study, we will stick to tuning two parameters, namely the mtry and the ntree parameters that have the following affect on our random forest model. 
#There are many other parameters, but these two parameters are perhaps the most likely to have the biggest effect on your final accuracy.

# we still use fitControl from previous analysis
#Direct from the help page for the randomForest() function in R:

#    mtry: Number of variables randomly sampled as candidates at each split.
#    ntree: Number of trees to grow.

#Let’s create a baseline for comparison by using the recommend defaults for each parameter and mtry=floor(sqrt(ncol(x))) or mtry around 19 and ntree=500.

mtry <- floor(sqrt(ncol(dat1.train)))
tunegrid <- expand.grid(.mtry=mtry)
rf1_default <- train(PHENOTYPE~., data=dat1.train, method="rf", metric="ROC", tuneGrid=tunegrid, trControl=fitControl)
test_pred_RF1_default <- predict(rf1_default, newdata = dat1.test)
confusionMatrix(test_pred_RF1_default,as.factor(dat1.test$PHENOTYPE),positive = "YES")
       
#Lets define a grid of algorithm parameters to try
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid",classProbs = TRUE)
set.seed(3233)
tunegrid <- expand.grid(.mtry=c(2:25))
rf1_gridsearch <- train(PHENOTYPE~., data=dat1.train, method="rf", metric="ROC", tuneGrid=tunegrid, trControl=control)
#print(rf1_gridsearch)
#plot(rf1_gridsearch)
#get the best mtry for the forest. 
best_mtry = rf1_gridsearch$finalModel$mtry
#Lets rerun with best_mtry but this time increase the forest to 1001 trees.( A suggested number of trees in the fiels is at least 1000)
set.seed(3233)
rf1_best_mtry_1001 <- train(PHENOTYPE~., data=dat1.train, method="rf", metric="ROC",ntree=1001, trControl=control,tuneGrid=data.frame(mtry=best_mtry))
print(rf1_best_mtry_1001)

test_pred_RF1_best_mtry_1001 <- predict(rf1_best_mtry_1001, newdata = dat1.test)
confusionMatrix(test_pred_RF1_best_mtry_1001,as.factor(dat1.test$PHENOTYPE),positive = "YES")

