# phenotype-prediction
R code to predict a phenotype based on genotypes, using SVM and Random Forest classification algorithms

This repository shows part of a project to evaluate the importance of genetic predisposition for a given clinical outcome. The dataset comes from a GWAS experiment where about 400 subjects were genotyped using an Illumina Infinium Global Screening array.

The first part of the project involves analyzing the GWAs data using Plink Software, following standard procedures. In addition to this, it also involves an imputation step with Beagle and further filtering afterwards. 

From here the full dataset is split in 5 folds for cross validation using scikit-learn and stratification using the phenotype proportion. For each of the cross validation folds, an additive model of inheritance is used to select associated variants. So we end up with five sets of associated genetic variants. Then the dataset is exported from Plink in text format to have five training sets and five validatoin sets with the associated variants at each fold. these datasets are imported into R for further analysis. These analysis are represented in this repository
