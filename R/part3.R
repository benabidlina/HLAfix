if(getRversion() >= "2.15.1") utils::globalVariables(".")
globalVariables(c("pheno", "Disease","P_allelic","P_dominant","Freq","freq","pdf","freq_control","freq_case","dev.off","AA","threshold","occurency", "outputFile", "HLA_KIR","capture.output","haploReg","Allele","val1","val2","Nctrl","Ncase","subdata","prob","Mean_pheno","ggplot2"))


#####################################################################
### Allele association analysis with trait(s) (-/+ Covariates) ######
#####################################################################

#' Analyze association of HLA alleles with many variate phenotypes (like case/control or quantitative trait)
#'
#' Take as input :
#' (1) a file containing the following 14 colunms :  subjectID Disease PC1 PC2  A A B B C C DRB1 DRB1 DQB1 DQB1 ("," separator between fields)
#' (2 optional) file with differents quantitative phenotypes : subjectID traitName1 traitName2 ... traitName(n)
#' (3 optional) file with covariable to correct regressions
#'
#' @param dataPath HLA alleles file
#' @param phenoPath phenotypes file with different trait's names in columns
#' @param covarPath must contains SubjectID followed by covariables
#' @param outputFolder folder where to save the result
#' @param outputFile result's filename
#' @param type either "disease" for case/control analysis or "pheno" with quantitative phenotype(s)
#' @return Linear Regression result and many plots : QQplot, personnalized manhattan plot and/or variation of phenotype's expression
#' @import ggplot2
#' @import plotrix
#' @import graphics
#' @importFrom  qqman qq
#' @importFrom stats binomial glm na.omit p.adjust lm
#' @importFrom utils read.table write.table capture.output
#' @importFrom stats binomial glm na.omit
#' @importFrom grDevices dev.off  pdf
#' @export
#'

hlaAlleleAnalysis <- function( dataPath , phenoPath=NULL , covarPath=NULL  , outputFolder , outputFile  , type = "pheno" )
{
  # Check Parameters
  if(!file.exists(dataPath)){stop("dataPath is not found ")}
  if(!dir.exists(outputFolder)){stop("outputFolder directoy doesn't exist. ")}
  if(!is.character(outputFile)){stop("outputFile must be character")}
  if( type=="pheno")
  {
    if(is.null(phenoPath)){stop("For type=pheno you must provide phenoPath parameter.")}
    if(!file.exists(phenoPath)){stop("phenoPath is not found ")}
  }

  ################  DATA IMPORTATION   ################
  all_data <-  as.data.frame(utils::read.table(dataPath  , stringsAsFactors = F, header = T, sep=','))
  if(ncol(all_data)!= 14){stop("data must contains 14 colunms :subjectID Disease PC1 PC2  A A B B C C DRB1 DRB1 DQB1 DQB1 ("," separator between fields)")}
  names(all_data) <- c("subjectID" , "Disease","PC1","PC2" ,"A","A.1","B","B.1","C","C.1","DRB1","DRB1.1","DQB1","DQB1.1")
  # if((!is.element(0,all_data$Disease) & !is.element(1,all_data$Disease))
  # {
  #   all_data$Disease <- ifelse(all_data$Disease=="Ctrl" || all_data$Disease=="ctrl" ||all_data$Disease=="control" || all_data$Disease=="Control" || all_data$Disease=="controle" ,0,1)  # create pheno and 0 for ctrl or 1 for cases
  #   if (sum(is.na(all_data$Disease)) == length(all_data$Disease)) {stop("Please code 0 and 1 the Disease statut of your individuals.")}
  # }
  #
  ################ OPTIONAL PARAMETERS  ################
  if(!is.null(covarPath)){
    covarCorrection <- stats::na.omit(utils::read.table(covarPath ,  header = TRUE, sep = "," , dec = "."))
    colnames(covarCorrection)[1] <- "subjectID"
    covarNames <- stats::na.omit(names(covarCorrection)) # extract name of covariates.
    covarNames <- covarNames[covarNames!="subjectID"]
    covariables <- noquote(paste(covarNames, collapse = " + "))
    all_data <- merge(all_data, covarCorrection, by ="subjectID")
    names(all_data) <- c("subjectID" , "Disease","PC1","PC2" ,"A","A.1","B","B.1","C","C.1","DRB1","DRB1.1","DQB1","DQB1.1", "AGE.V0", "SEX","CMV.V0","NBYTABAC")
  }

  # THRESHOLDS
  freq_threshold <- readline(prompt = "\n \n Enter a frequency threshold of alleles to displays in manhattan plot\n Frequency : " )
  threshold <- readline(prompt = "In order to display less alleles and more clearly. Define a p value treshold under which one alleles are not displayed \n P-value:  ")

  switch(type,
         disease={

           HLA_genes <- c("A", "B", "C", "DRB1", "DQB1") # Set up HLA genes
           res <-  c("Locus","HLA_allele", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant") ## Create header to store results

           # File that will contains best Pvalue for all phenotypes
           X <- c("Locus","HLA_allele", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")
           write.table(t(X), file=paste(outputFolder,"BestPvalue_AllelesDiseaseAnalysis.csv ",sep = "") , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")

           for(gene in HLA_genes) # Make a loop with all the HLA genes ## Create the loop to analyze the data one gene at a time and inside another loop that takes each alleles
           {
             subdata <- na.omit(all_data[,c("subjectID","Disease","PC1","PC2",gene,paste(gene,".1",sep = ""))])
             names(subdata) <- c("subjectID","Disease","PC1","PC2", "allele1" , "allele2" )
             if(is.null(subdata) | length(subdata) == 0) {stop(" \n No rows raimining in your dataset, too much NA. Is your disease statut coded 0 and 1 ? \n ")}
             HLA <- stats::na.omit(as.character(unique(c(unique(subdata$allele1), unique(subdata$allele2))))) # Make a list of all existing HLA alleles in each file
             for(allele in HLA)
             { # Make a loop with all the HLA alleles for each gene, Convert the allele to 1,2 code (val1, val2), then to genotype (val3), and finally to recessive (val4) or dominant (val5) testable model
               val1 <- paste("allele1_", allele, sep="")
               subdata[,val1] <- ifelse(subdata$allele1 == allele, 1,
                                        ifelse(subdata$allele1 != allele, 2, NA))

               val2 <- paste("allele2_", allele, sep="")
               subdata[,val2] <- ifelse(subdata$allele2 == allele, 1,
                                        ifelse(subdata$allele2 != allele, 2, NA))

               val3 <- paste("genotype_", allele, sep="")
               subdata[,val3] <- paste(subdata[,val1], subdata[,val2], sep="")

               val4 <- paste("genoA_", allele, sep="")
               subdata[,val4] <- ifelse(subdata[,val3] == '11', 1,
                                        ifelse(subdata[,val3] == '12', 0,
                                               ifelse(subdata[,val3] == '21', 0,
                                                      ifelse(subdata[,val3] == '22', -1, NA))))

               val5 <- paste("genoD_", allele, sep="")
               subdata[,val5] <- ifelse(subdata[,val3] == '11', 1,
                                        ifelse(subdata[,val3] == '12', 1,
                                               ifelse(subdata[,val3] == '21', 1,
                                                      ifelse(subdata[,val3] == '22', 0, NA))))


               freq <- sum(subdata[,val1]==1,subdata[,val2]==1)/(nrow(subdata)*2) # Overall frequency
               n_ind_ctrl <- nrow(subset(subdata, Disease == 0))
               n_ind_cas <-nrow(subset(subdata, Disease == 1))

               if( n_ind_ctrl != 0 & n_ind_cas != 0 ) # Frequency are calculated for case and control
               {
                 # N <- 4/((1/nrow(subset(subdata, Disease == 1)))+(1/nrow(subset(subdata, Disease == 0))))
                 N <- nrow(subdata[,!is.na(Disease)])
                 freq_control <- sum(subset(subdata, Disease == 0)[,val1]==1,subset(subdata, Disease == 0)[,val2]==1)/(nrow(subset(subdata, Disease == 0))*2) # Controls frequency
                 freq_case <- sum(subset(subdata, Disease == 1)[,val1]==1,subset(subdata, Disease == 1)[,val2]==1)/(nrow(subset(subdata, Disease == 1))*2) # Cases frequency with homozygote for the allele (11)

               }else{
                 if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl and set to 0
                   freq_control <- 0
                   freq_case <- sum(subset(subdata, Disease == 1)[,val1]==1,subset(subdata, Disease == 1)[,val2]==1)/(nrow(subset(subdata, Disease == 1))*2) # Cases frequency with homozygote for the allele (11)
                   N <- nrow(subset(subdata, Disease == 1))

                 }else{
                   if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for ctrl and set to 0
                     freq_control <- sum(subset(subdata, Disease == 0)[,val1]==1,subset(subdata, Disease == 0)[,val2]==1)/(nrow(subset(subdata, Disease == 0))*2) # Controls frequency
                     freq_case <- 0
                     N <- nrow(subset(subdata, Disease == 0))

                   }else{ # Unknown case control Status
                     stop("For a Case/Control analysis, the Disease colunm is required. Please note 1=case and 0=control")

                   }
                 }
               }

               # Adapt the regression correction co-variables
               subdata[, val4] <- as.factor(subdata[, val4]) # convert genoA and genoD to factor for a logistic regression
               subdata[, val5] <- as.factor(subdata[, val5])
               subdata[,"Disease"] <- as.factor(subdata[,"Disease"])

               if(!is.null(covarPath)){
                 reg_allelic <- summary(stats::glm(Disease ~ subdata[, val4] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = subdata, family=stats::binomial(link="logit"),maxit=100)) # Logistic regression with allelic model
                 reg_dominant <- summary(stats::glm(Disease ~ subdata[, val5] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = subdata, family=stats::binomial(link="logit"),maxit=100)) # Logistic regression with dominant model

               }else{
                 reg_allelic <- summary(stats::glm(Disease ~ subdata[, val4] + PC1 + PC2, data = subdata, family=stats::binomial(link="logit"),maxit=100)) # Logistic regression with allelic model
                 reg_dominant <- summary(stats::glm(Disease ~ subdata[, val5] + PC1 + PC2, data = subdata, family=stats::binomial(link="logit"),maxit=100)) # Logistic regression with dominant model

               }

               N <- 4/((1/nrow(subset(subdata, Disease == 1)))+(1/nrow(subset(subdata, Disease == 0))))
               P_allelic <- reg_allelic$coefficients[2,4]
               OR_allelic <- exp(reg_allelic$coefficients[2,1])
               lower_95_CI_allelic <- exp(reg_allelic$coefficients[2,1]-1.96*reg_allelic$coefficients[2,2])
               upper_95_CI_allelic <- exp(reg_allelic$coefficients[2,1]+1.96*reg_allelic$coefficients[2,2])
               Beta_allelic <- reg_allelic$coefficients[2,1]
               SE_allelic <- reg_allelic$coefficients[2,2]
               P_dominant <- reg_dominant$coefficients[2,4]
               OR_dominant <- exp(reg_dominant$coefficients[2,1])
               lower_95_CI_dominant <- exp(reg_dominant$coefficients[2,1]-1.96*reg_dominant$coefficients[2,2])
               upper_95_CI_dominant <- exp(reg_dominant$coefficients[2,1]+1.96*reg_dominant$coefficients[2,2])
               Beta_dominant <- reg_dominant$coefficients[2,1]
               SE_dominant <- reg_dominant$coefficients[2,2]

               # add the regression result to the "res" vector
               res <- rbind(res, c(gene, allele , freq , freq_control , freq_case, N, P_allelic, OR_allelic, lower_95_CI_allelic, upper_95_CI_allelic, Beta_allelic, SE_allelic, P_dominant, OR_dominant, lower_95_CI_dominant, upper_95_CI_dominant, Beta_dominant, SE_dominant)) # See header of res for elements explanation

             }

           }
           # Save results
           # if(!is.null(res) & length(res) > 18){
           utils::write.table(res, paste(outputFolder , outputFile,"DiseaseAnalysis.csv" , sep=""), row.names=F, quote=F, col.names=F, sep = ",")
           cat("Regression at allele levels is processed successfully \n")
           cat(paste("File", paste(outputFolder , outputFile,"DiseaseAnalysis.csv" , sep="") ,'is created.  \n', sep = " "))

           ##########################################################################################################################################
           information <-   utils::read.table( paste(outputFolder , outputFile,"DiseaseAnalysis.csv" , sep=""), stringsAsFactors = F, header = T, sep=',')
           seuil <- pvalueCorrection(alleleRegResult = paste(outputFolder , outputFile,"DiseaseAnalysis.csv" , sep=""), correction = "bonferroni")
           if(is.null(seuil)){stop("Seuil calculated with pvalueCorrection() is null")}

           grDevices::pdf(file = paste(outputFolder ,outputFile,"Alleles_Disease_Analysis_freq",freq_threshold,"_pvalManhattan",threshold , ".pdf" ,sep=""))
           ################### Best Pvalue selection gathered in an unique file
           resultLoop <- as.data.frame(information)
           resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold  & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05) ,]
           if(!is.null(resultLoop) & length(resultLoop)> 0){
             suppressWarnings({
               resultLoop <-resultLoop[,c("Locus","HLA_allele", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")]
               write.table(resultLoop, file=paste(outputFolder,"BestPvalue_AllelesDiseaseAnalysis.csv ",sep = "") , append = TRUE, row.names = FALSE ,  col.names = FALSE ,quote=F, sep = ",")
             })
           }
           ################### QQPLOT
           ## Lambda: genomic inflation factor, show the presence of bias
           # allelic Model
           observedR <- sort(information$P_allelic)
           lobsR <- -(log10(observedR))
           expectedR <- c(1:length(observedR))
           lexpR <- -(log10(expectedR / (length(expectedR)+1)))
           lambdaR<- summary(stats::lm(lobsR~lexpR))

           # Dominant Model
           observedD <- sort(information$P_dominant)
           lobsD <- -(log10(observedD))
           expectedD <- c(1:length(observedD))
           lexpD <- -(log10(expectedD / (length(expectedD)+1)))
           lambdaD<- summary(stats::lm(lobsD~lexpD))

           # QQ PLOT
           qq1 <- qq(information$P_allelic , main="Disease analysis of alleles in allelic model ", sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
           qq2 <- qq(information$P_dominant , main="Disease analysis of alleles in dominant model ", sub= paste("Lambda = ", round(lambdaD$coefficients[2,1],2),sep = ""))
           capture.output( print(qq1),  file='NUL' )
           capture.output( print(qq2) , file='NUL')

           ################### MANHATTAN PLOT
           selection <-  na.omit(information[ information$Frequency >= freq_threshold & (information$P_allelic<=threshold | information$P_dominant<=threshold),c(1,2,7,13)]) # "Locus" "HLA_allele" "P_allelic" "P_dominant"
           if(nrow(selection)==0){
             cat(paste("None alleles to display. You can change your threshold p-value. Set at ",threshold,sep = ""))
           }else{
             options(scipen = 999)
             colnames(selection) <- c("Locus", "Allele" , "P_allelic", "P_dominant")
             maxYvalPallelic <- (max(-log10(selection$P_allelic))+1)
             maxYvalPdom <- (max(-log10(selection$P_dominant))+1)

             # Adapt x axis legend with the number of alleles to plot
             n <- nrow(selection)
             if(n<50){
               size_txt_x_axis <- 10
               size_point <- 2
             }else{
               if(n>50 & n<95){
                 size_txt_x_axis <- 5.5
                 size_point <- 2
               }else{
                 if(n>95){
                   size_txt_x_axis <- 3.5
                   size_point <- 1
                 }
               }
             }
             plot1 <-ggplot( selection , aes(x=Allele, y=  -(log10(P_allelic)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(selection$P_allelic)) >= seuil, "red","black")) + scale_y_continuous( limits = c(0 , maxYvalPallelic)) +  theme(axis.text.x = element_text(size = size_txt_x_axis, angle=90)) + facet_grid(.~Locus , scales = "free" , space = "free_x") + geom_hline(yintercept= -log10(0.05) , color = "blue")  + geom_hline(yintercept= seuil , color = "red")  + ggtitle(paste("Disease analysis of alleles in allelic model",sep = ""))
             plot2 <-ggplot( selection , aes(x=Allele, y=  -(log10(P_dominant)) ))  + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(selection$P_dominant)) >= seuil, "red","black")) + scale_y_continuous( limits = c(0 , maxYvalPdom)) +  theme(axis.text.x = element_text(size = size_txt_x_axis, angle=90)) + facet_grid(.~Locus , scales = "free" , space = "free_x") + geom_hline(yintercept= -log10(0.05) , color = "blue")  + geom_hline(yintercept= seuil , color = "red") + ggtitle(paste("Disease analysis of alleles in dominant model",sep = ""))
           }
           suppressWarnings({
             print(plot1)
             print(plot2)
           })

           grDevices::dev.off()

         }, # end case = disease
  pheno ={
           #List of phenotypes
           phenoFile <-  as.data.frame(utils::read.table(phenoPath  , stringsAsFactors = F, header = T, sep=','))
           colnames(phenoFile)[1] <- "subjectID"
           phenotypes <- stats::na.omit(unique(c(names(phenoFile)))) # Make a list of of phenotypes in the file from colnames
           phenotypes <- phenotypes[phenotypes != "subjectID"]

           # File that will contains best Pvalue for all phenotypes
           X <- c("Phenotype","Locus","HLA_allele", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")
           write.table(t(X), file=paste(outputFolder,"BestPvalue_AllelesPhenoAnalysis.csv ",sep = "") , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")
           grDevices::pdf(file = paste(outputFolder ,outputFile,"Alleles_Pheno_Analysis_freq",freq_threshold,"_pvalManhattan",threshold , ".pdf" ,sep=""))

           for(pheno in phenotypes)
           {
             print(paste(pheno, " is currently being analysed.", sep = ""))
             HLA_genes <- c("A", "B", "C", "DRB1", "DQB1")
             res <- c("Locus","HLA_allele", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")
             phenoGeno <- merge(all_data,phenoFile[,c("subjectID",pheno)],by = "subjectID") # could contain co-variables if provided
             colnames(phenoGeno)[ncol(phenoGeno)] <- "trait"
             phenoGeno <- na.omit(phenoGeno)
             Mean_pheno <- mean(phenoGeno$trait)

             if(nrow(phenoGeno)==0){stop("After merging genetic data with pheno and removing NA there is no data left to run an analysis. Maybe an incompatible sample.ID")}
             if(dir.exists(paste(outputFolder,"expression_HLA/",sep = ""))){
               grDevices::pdf(file = paste(outputFolder ,"expression_HLA/",outputFile,pheno,"-expression-vs-alleles " , ".pdf" ,sep=""))
             }else{
               dir.create(paste(outputFolder,"expression_HLA/",sep = ""))
               grDevices::pdf(file = paste(outputFolder ,"expression_HLA/",outputFile,pheno,"-expression-vs-alleles " , ".pdf" ,sep=""))
             }
             for(gene in HLA_genes) # Make a loop with all the HLA genes
             {
               if(!is.null(covarPath)){
                 subdata <- na.omit (subset(phenoGeno, select=c("subjectID","Disease","PC1","PC2",gene,paste(gene,".1",sep = ""),"trait", "AGE.V0", "SEX","CMV.V0","NBYTABAC")))
                 names(subdata) <- c("subjectID","Disease", "PC1","PC2","allele1" , "allele2", "trait" , "AGE.V0", "SEX","CMV.V0","NBYTABAC")
               }else{
                 subdata <- na.omit (subset(phenoGeno, select=c("subjectID","Disease","PC1","PC2",gene,paste(gene,".1",sep = ""),"trait")))
                 names(subdata) <- c("subjectID","Disease", "PC1","PC2","allele1" , "allele2", "trait")
               }
               HLA <- stats::na.omit(as.character(unique(c(unique(subdata$allele1), unique(subdata$allele2))))) # Make a list of all existing HLA alleles in each file
               for(allele in HLA)
               {
                 val1 <- paste("allele1_", allele, sep="")
                 subdata[,val1] <- ifelse(subdata$allele1 == allele, 1,
                                          ifelse(subdata$allele1 != allele, 2, NA))

                 val2 <- paste("allele2_", allele, sep="")
                 subdata[,val2] <- ifelse(subdata$allele2 == allele, 1,
                                          ifelse(subdata$allele2 != allele, 2, NA))

                 val3 <- paste("genotype_", allele, sep="")
                 subdata[,val3] <- paste(subdata[,val1], subdata[,val2], sep="")

                 val4 <- paste("genoA_", allele, sep="")
                 subdata[,val4] <- ifelse(subdata[,val3] == '11', 1,
                                          ifelse(subdata[,val3] == '12', 0,
                                                 ifelse(subdata[,val3] == '21', 0,
                                                        ifelse(subdata[,val3] == '22', -1, NA))))

                 val5 <- paste("genoD_", allele, sep="")
                 subdata[,val5] <- ifelse(subdata[,val3] == '11', 1,
                                          ifelse(subdata[,val3] == '12', 1,
                                                 ifelse(subdata[,val3] == '21', 1,
                                                        ifelse(subdata[,val3] == '22', 0, NA))))


                 freq <- sum(subdata[,val1]==1,subdata[,val2]==1)/(nrow(subdata)*2) # Overall frequency
                 n_ind_ctrl <- nrow(subset(subdata, Disease == 0))
                 n_ind_cas <-nrow(subset(subdata, Disease == 1))

                 if( n_ind_ctrl != 0 & n_ind_cas != 0 ) # Frequency are calculated for case and control
                 {
                   # N <- 4/((1/nrow(subset(subdata, Disease == 1)))+(1/nrow(subset(subdata, Disease == 0))))
                   N <- nrow(subdata[,!is.na(Disease)])
                   freq_control <- sum(subset(subdata, Disease == 0)[,val1]==1,subset(subdata, Disease == 0)[,val2]==1)/(nrow(subset(subdata, Disease == 0))*2) # Controls frequency
                   freq_case <- sum(subset(subdata, Disease == 1)[,val1]==1,subset(subdata, Disease == 1)[,val2]==1)/(nrow(subset(subdata, Disease == 1))*2) # Cases frequency with homozygote for the allele (11)

                 }else{
                   if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl and set to 0
                     freq_control <- 0
                     freq_case <- sum(subset(subdata, Disease == 1)[,val1]==1,subset(subdata, Disease == 1)[,val2]==1)/(nrow(subset(subdata, Disease == 1))*2) # Cases frequency with homozygote for the allele (11)
                     N <- nrow(subset(subdata, Disease == 1))

                   }else{
                     if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for ctrl and set to 0
                       freq_control <- sum(subset(subdata, Disease == 0)[,val1]==1,subset(subdata, Disease == 0)[,val2]==1)/(nrow(subset(subdata, Disease == 0))*2) # Controls frequency
                       freq_case <- 0
                       N <- nrow(subset(subdata, Disease == 0))

                     }else{ # Unknown case control Status
                       N <- 0
                       freq_control <- 0
                       freq_case <- 0
                     }
                   }
                 }

                 # Adapt the regression correction co-variables
                 if(!is.null(covarPath)){
                   reg_allelic <- summary(stats::lm(trait ~ subdata[, val4] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = subdata  )) # Linear regression with allelic model
                   reg_dominant <- summary(stats::lm(trait ~ subdata[, val5] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = subdata )) # Linear regression with dominant model

                 }else{
                   reg_allelic <- summary(stats::lm(trait ~ subdata[, val4] + PC1 + PC2, data = subdata)) # Linear regression with allelic model
                   reg_dominant <- summary(stats::lm(trait ~ subdata[, val5] + PC1 + PC2, data = subdata)) # Linear regression with dominant model

                 }

                 P_allelic <- reg_allelic$coefficients[2,4]
                 Beta_allelic <- round(reg_allelic$coefficients[2,1],2)
                 SE_allelic <- round(reg_allelic$coefficients[2,2],2)
                 Rsquare_allelic <-round(reg_allelic$adj.r.squared, 3)

                 P_dominant <- reg_dominant$coefficients[2,4]
                 Beta_dominant <- round(reg_dominant$coefficients[2,1],2)
                 SE_dominant <- round(reg_dominant$coefficients[2,2],2)
                 Rsquare_dominant <-round(reg_dominant$adj.r.squared, 3)
                 # add the regression result to the "res" vector
                 options(scipen = 999)
                 res <- rbind(res, c(gene, allele , freq, freq_control, freq_case, N, Mean_pheno, P_allelic, Beta_allelic, SE_allelic, Rsquare_allelic, P_dominant, Beta_dominant, SE_dominant, Rsquare_dominant)) # See header of res for elements explanation

               }

               ################  ehPlot for visualization of quantitative phenotype as a function of allele ################
               selection <- phenoGeno[,c("subjectID" , gene, paste(gene,".1",sep=""),"trait")]
               names(selection) <- c("subjectID" , "allele1","allele2","trait")
               selection$genotype <- ifelse(selection$allele1==selection$allele2,2,1) # coded 1:heterozugous 2:homozygous
               selection$genotype <- as.factor(selection$genotype)
               names(selection) <- c("subjectID" , "allele1","allele2","trait","genotype")

               # alleles independally
               HLA_a1 <- na.omit(selection[,c("trait","allele1", "genotype")]) #allele 1 without NA
               names(HLA_a1)[2] <- "HLA_alleles" #rename column
               HLA_a2 <- na.omit(selection[,c("trait","allele2", "genotype")]) #allele 2 without NA
               names(HLA_a2)[2] <- "HLA_alleles" #rename column
               HLA_all_alleles <- rbind(HLA_a1, HLA_a2) #combine data
               names(HLA_all_alleles) <- c("trait","All","genotype")
               HLA_all_alleles$HLA_alleles_2d <- substr(HLA_all_alleles$All, 1, 2) # Create 2 digits column

               #Remove alleles with less than 5 data point
               alleles <- unique(HLA_all_alleles$HLA_alleles_2d)
               for(i in alleles){
                 if(sum(HLA_all_alleles$HLA_alleles_2d == i) < 5){
                   HLA_all_alleles <- HLA_all_alleles[! HLA_all_alleles$HLA_alleles_2d == i,]
                 }
               }
               if(nrow(HLA_all_alleles)!=0)
               { # try to draw the plot only if data are available
                 #Get the stat
                 HLA_all_alleles_boxplot <- boxplot(HLA_all_alleles$trait ~ HLA_all_alleles$HLA_alleles_2d, plot=F)
                 HLAstat <- HLA_all_alleles_boxplot$stats
                 colnames(HLAstat) <- HLA_all_alleles_boxplot$names
                 rownames(HLAstat) <- c('min','lower quartile','median','upper quartile','max')
                 HLAstat <- t(HLAstat)
                 HLAstat <- HLAstat[order(HLAstat[,3]),]
                 HLAstat <- cbind(HLAstat, as.matrix(tapply(HLA_all_alleles$trait, HLA_all_alleles$HLA_alleles_2d, mean)))
                 colnames(HLAstat)[6] <- "mean"
                 # write.table(HLAstat, "HLA_DR_2digits_expr_MI_data_280318.txt", sep="\t")
                 HLA_all_alleles$HLA_alleles_2d <- factor(HLA_all_alleles$HLA_alleles_2d, levels=rownames(HLAstat)) #set the HLA-C alleles as factors to show the graph in increasing order

                 #Draw the graph
                 suppressWarnings({
                   myEhplot <- plotrix::ehplot(HLA_all_alleles$trait, HLA_all_alleles$HLA_alleles_2d, xlab=paste("HLA-", gene," alleles",sep = ""), ylab=paste(pheno," expression"), main=paste(pheno," expression depending on alleles"), pch=20, interval=5, col=HLA_all_alleles$genotype)
                   capture.output( print(myEhplot) , file='NUL')
                 })
               }
             } # end gene loop
             grDevices::dev.off()

             # save a result file for one phenotype.
             utils::write.table(res, paste(outputFolder , outputFile , "-", pheno , ".csv" ,sep=""), row.names=F, quote=F, col.names=F, sep = ",")
             cat(paste("Regression at allele levels and with : ", pheno , " as covariates  is processed successfully \n" , sep = ""))
             cat( paste("File", paste(outputFolder , outputFile, "-", pheno   , sep="") ,'is created.  \n', sep = " "))

             ##########################################################################################################################################

             information <-   utils::read.table(paste(outputFolder , outputFile , "-", pheno , ".csv" ,sep=""), stringsAsFactors = F, header = T, sep=',')
             seuil <- pvalueCorrection(alleleRegResult = paste(outputFolder , outputFile , "-", pheno , ".csv" ,sep=""), correction = "bonferroni")
             if(is.null(seuil)){stop("Seuil calculated with pvalueCorrection() is null")}

             # grDevices::dev.off()
             ################### Best Pvalue selection gathered in an unique file
             resultLoop <- as.data.frame(information)
             resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold  & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05) ,]
             if(nrow(resultLoop)!=0){
               suppressWarnings({
                 options(scipen = 999)
                 resultLoop$Phenotype <- pheno
                 resultLoop <-resultLoop[,c("Phenotype","Locus","HLA_allele", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")]
                 write.table(resultLoop, file=paste(outputFolder,"BestPvalue_AllelesPhenoAnalysis.csv ",sep = "") , append = TRUE, row.names = FALSE ,  col.names = FALSE ,quote=F, sep = ",")
               })
             }
             ################### QQPLOT

             ## Lambda: genomic inflation factor, show the presence of bias
             # allelic Model
             observedR <- sort(information$P_allelic)
             lobsR <- -(log10(observedR))
             expectedR <- c(1:length(observedR))
             lexpR <- -(log10(expectedR / (length(expectedR)+1)))
             lambdaR<- summary(stats::lm(lobsR~lexpR))

             # Dominant Model
             observedD <- sort(information$P_dominant)
             lobsD <- -(log10(observedD))
             expectedD <- c(1:length(observedD))
             lexpD <- -(log10(expectedD / (length(expectedD)+1)))
             lambdaD<- summary(stats::lm(lobsD~lexpD))

             # QQ PLOT
             qq1 <- qq(information$P_allelic , main=paste(pheno," in allelic modele ",sep = ""), sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
             qq2 <- qq(information$P_dominant , main=paste(pheno,"Alleles in Dominant modele ",sep = ""), sub= paste("Lambda = ", round(lambdaD$coefficients[2,1],2),sep = ""))
             capture.output( print(qq1),  file='NUL' )
             capture.output( print(qq2) , file='NUL')

             ################### MANHATTAN PLOT
             selection <-  information[ information$Frequency >= freq_threshold & (information$P_allelic <= threshold | information$P_dominant <= threshold),c(1,2,8,12)] # "Locus" "HLA_allele" "P_allelic" "P_dominant"
             if(nrow(selection)==0){
               cat(paste("None alleles to display. You can change your threshold p-value. Set at ",threshold,sep = ""))
             }else{
               options(scipen = 999)
               colnames(selection) <- c("Locus", "Allele" , "P_allelic", "P_dominant")
               maxYvalPallelic <- (max(-log10(selection$P_allelic))+1)
               maxYvalPdom <- (max(-log10(selection$P_dominant))+1)

               # Adapt x axis legend with the number of alleles to plot
               n <- nrow(selection)
               if(n<50){
                 size_txt_x_axis <- 10
                 size_point <- 2
               }else{
                 if(n>50 & n<95){
                   size_txt_x_axis <- 5.5
                   size_point <- 2
                 }else{
                   if(n>95){
                     size_txt_x_axis <- 3.5
                     size_point <- 1
                   }
                 }
               }
               plot1 <-ggplot( selection , aes(x=Allele, y=  -(log10(P_allelic)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(selection$P_allelic)) > seuil, "red","black")) + scale_y_continuous( limits = c(0 , maxYvalPallelic)) +  theme(axis.text.x = element_text(size = size_txt_x_axis, angle=90)) + facet_grid(.~Locus , scales = "free" , space = "free_x") + geom_hline(yintercept= -log10(0.05) , color = "blue")  + geom_hline(yintercept= seuil , color = "red")  + ggtitle(paste( pheno ," analysis in allelic Model",sep = ""))
               plot2 <-ggplot( selection , aes(x=Allele, y=  -(log10(P_dominant)) ))  + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(selection$P_dominant)) > seuil, "red","black")) + scale_y_continuous( limits = c(0 , maxYvalPdom)) +  theme(axis.text.x = element_text(size = size_txt_x_axis, angle=90)) + facet_grid(.~Locus , scales = "free" , space = "free_x") + geom_hline(yintercept= -log10(0.05) , color = "blue")  + geom_hline(yintercept= seuil , color = "red") + ggtitle(paste( pheno ," analysis in Dominant Model",sep = ""))
             }
             suppressWarnings({
               print(plot1)
               print(plot2)
             })
           } # end pheno loop

           grDevices::dev.off()
         }) # end case = pheno
}


#####################################################################
### Haplotype association analysis with trait(s) (-/+ Covariates)###
#####################################################################

#' Analyze association of HLA haplotypes with many variate phenotypes (like case/control or quantitative trait)
#' Take as input:
#' (1 required) haplotypes :  subjectID haplotype1 freq1  haplotype2 freq2 postprob ("," separator between fields)
#' (2 optional) differents quantitative phenotypes : subjectID traitName1 traitName2 ... traitName(n)
#' (3 required) Disease information : 0= control 1=case (needed for frequencies and/or type=disease analysis)
#' (4 required) the 2 first principals components of the PCA
#' (5 optional) covariable to correct regressions
#'
#'
#' @param haploDataPath HLA haplotypes data exactly as the one provided by easy-HLA webtool
#' @param phenoPath phenotypes file with different trait's names in columns
#' @param statut file conaining subject ID and Disease statut
#' @param pcaPath pca's result: subjectID, PC1 and PC2
#' @param covarPath must contains SubjectID followed by covariables
#' @param outputFolder folder where to save the result
#' @param outputFile result's filename
#' @param type  either "disease" for case/control analysis or "pheno" with quantitative phenotype(s)
#' @return Linear Regression result and many plots : QQplot and personnalized manhattan plot
#' @import ggplot2
#' @import ggrepel
#' @importFrom  qqman qq
#' @importFrom stats binomial glm na.omit p.adjust lm
#' @importFrom utils read.table write.table capture.output
#' @importFrom grDevices dev.off pdf
#' @export
#'


hlaHaplotypeAnalysis <- function( haploDataPath , statut , phenoPath= NULL, pcaPath=NULL, covarPath=NULL , outputFolder , outputFile, type="pheno")
{
  # Check Parameters
  if(!file.exists(haploDataPath)){stop("haploDataPath is not found ")}
  if(!file.exists(phenoPath)){stop("phenoPath is not found ")}
  if(!file.exists(statut)){stop("statut is not found ")}
  if(!dir.exists(outputFolder)){stop("outputFolder directoy doesn't exist. ")}
  if(!is.character(outputFile)){stop("outputFile must be character")}

  if( type=="pheno")
  {
    if(is.null(phenoPath)){stop("For type=pheno you must provide phenoPath parameter.")}
    if(!file.exists(phenoPath)){stop("phenoPath is not found ")}
  }

  ################  DATA IMPORTATION   ################
  haploData <- as.data.frame(utils::read.table(haploDataPath  , stringsAsFactors = F, header = T, sep=','))
  phenoFile <-  as.data.frame(utils::read.table(phenoPath  , stringsAsFactors = F, header = T, sep=','))
  statut <- stats::na.omit(utils::read.table(statut , header = TRUE , sep = "," , dec = "."))
  if(ncol(statut)!=2){stop("statut must contains 2 colunms : subjectID and Affection statut")}
  names(statut) <- c("subjectID","Disease")
  # statut$Disease <- ifelse(statut$Disease=="Ctrl" || statut$Disease=="ctrl" ||statut$Disease=="control" || statut$Disease=="Control" || statut$Disease=="controle" ,0,1) # create pheno and 0 for ctrl or 1 for cases

  names(haploData) <- c("subjectID" ,"A.B.C.DR.DQ.1","freq1","A.B.C.DR.DQ.2","freq2","postProb" ,"precision")
  colnames(phenoFile)[1] <- "subjectID"

  all_data <- merge(haploData, statut, by ="subjectID")
  colnames(all_data)[ncol(all_data)] <- "Disease"

  ################ OPTIONAL PARAMETERS  ################
  if(!is.null(pcaPath)){
    PCs <- stats::na.omit(utils::read.table(pcaPath ,  header = TRUE, sep = "," , dec = "."))
    if(ncol(PCs)!=3){stop("PCs must contains 3 colunms : subjectID, PC1 and PC2.")}
    colnames(PCs)[1]<-"subjectID"
    colnames(PCs)[2]<-"PC1"
    colnames(PCs)[3]<-"PC2"
    all_data <- merge(all_data, PCs[,1:3], by ="subjectID")
  }
  if(!is.null(covarPath)){
    covarCorrection <- stats::na.omit(utils::read.table(covarPath ,  header = TRUE, sep = "," , dec = "."))
    colnames(covarCorrection)[1] <- "subjectID"
    all_data <- na.omit(merge(all_data, covarCorrection ,by = "subjectID"))
    colnames(all_data)[ncol(all_data)-1] <- "NBYTABAC"
    colnames(all_data)[ncol(all_data)-2] <- "CMV.V0"
    colnames(all_data)[ncol(all_data)-3] <- "SEX"
    colnames(all_data)[ncol(all_data)-4] <- "AGE.V0"

    # covarNames <- stats::na.omit(names(covarCorrection))  # extract name of covariates.
    # covarNames <- covarNames[covarNames!="subjectID"]
    # covariables <- noquote(paste(covarNames, collapse = " + "))
  }

  # THRESHOLDS
  freq_threshold <- readline(prompt = "\n \n Enter a frequency threshold of haplotype to displays in manhattan plot\n Frequency: " )
  threshold <- readline(prompt = "\n \n In order to display most interesting haplotypes and in clearly way. \n Define a p value treshold under which one haplotypes are not displayed in manhattan plot \n P-value:  ")
  occurrency <- readline(prompt = "\n \n Enter number of haplotype occurency in the cohort \n Occurency: " )

  switch(type,
         disease={
           res <- c("HLA_haplotype", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_95%_CI_dominant", "upper_95%_CI_dominant", "Beta_dominant", "SE_dominant")

           ## Determine haplotypes to analyze; more than 2 occurences
           df <- as.data.frame(table(stats::na.omit(c( as.character(all_data$A.B.C.DR.DQ.1) , as.character(all_data$A.B.C.DR.DQ.2))))) ## number of occurences
           haplo <-subset(df , Freq >= occurency) ## select haplo that appears more than 2 times
           grDevices::pdf(file = paste(outputFolder ,outputFile,"Haplo_Disease_analysis_freq",freq_threshold,"_occurency", occurrency ,"_pvalManhattan",threshold , ".pdf" ,sep=""))

           for(hap in haplo[,1])
           { # Make a loop with all the HLA haplotypes, Convert the haplotypes to 1,2 code (val1, val2), then to genotype (val3), and finally to recessive (val4) or dominant (val5) testable model
             val1 <- paste("allele1_", hap, sep="")
             all_data[,val1] <- ifelse(all_data$A.B.C.DR.DQ.1 == hap, 1,
                                       ifelse(all_data$A.B.C.DR.DQ.1 != hap, 2, NA))

             val2 <- paste("allele2_", hap, sep="")
             all_data[,val2] <- ifelse(all_data$A.B.C.DR.DQ.2 == hap, 1,
                                       ifelse(all_data$A.B.C.DR.DQ.2 != hap, 2, NA))

             val3 <- paste("genotype_", hap, sep="")
             all_data[,val3] <- paste(all_data[,val1], all_data[,val2], sep="")

             val4 <- paste("genoA_", hap, sep="")
             all_data[,val4] <- ifelse(all_data[,val3] == '11', 1,
                                       ifelse(all_data[,val3] == '12', 0,
                                              ifelse(all_data[,val3] == '21', 0,
                                                     ifelse(all_data[,val3] == '22', -1, NA))))

             val5 <- paste("genoD_", hap, sep="")
             all_data[,val5] <- ifelse(all_data[,val3] == '11', 1,
                                       ifelse(all_data[,val3] == '12', 1,
                                              ifelse(all_data[,val3] == '21', 1,
                                                     ifelse(all_data[,val3] == '22', 0, NA))))


             ## Frequencies
             freq <- sum(all_data[,val1]==1,all_data[,val2]==1)/(nrow(all_data)*2) # Overall frequency
             n_ind_ctrl <- nrow(subset(all_data, Disease == 0))
             n_ind_cas <-nrow(subset(all_data, Disease == 1))

             if( n_ind_ctrl != 0 & n_ind_cas != 0 )
             { # Frequency are calculated for case and control
               # N <- 4/((1/nrow(subset(all_data, Disease == 1)))+(1/nrow(subset(all_data, Disease == 0))))
               N <- nrow(subdata[,!is.na(Disease)])
               freq_control <- sum(subset(all_data, Disease == 0)[,val1]==1,subset(all_data, Disease == 0)[,val2]==1)/(nrow(subset(all_data, Disease == 0))*2) # Controls frequency
               freq_case <- sum(subset(all_data, Disease == 1)[,val1]==1,subset(all_data, Disease == 1)[,val2]==1)/(nrow(subset(all_data, Disease == 1))*2) # Cases frequency

             }else{
               if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl and set to 0
                 freq_control <- 0
                 freq_case <- sum(subset(all_data, Disease == 1)[,val1]==1,subset(all_data, Disease == 1)[,val2]==1)/(nrow(subset(all_data, Disease == 1))*2) # Cases frequency
                 N <- nrow(subset(all_data, Disease == 1))

               }else{
                 if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for case and set to 0
                   freq_control <- sum(subset(all_data, Disease == 0)[,val1]==1,subset(all_data, Disease == 0)[,val2]==1)/(nrow(subset(all_data, Disease == 0))*2) # Controls frequency
                   freq_case <- 0
                   N <- nrow(subset(all_data, Disease == 0))
                 }else{ # Unknown case control Status
                   stop("For a Case/Control analysis, the Disease colunm is required. Please note 1=case and 0=control")
                 }
               }
             }
             # Adapt the regression correction co-variables
             if(is.null(covarPath) & is.null(pcaPath)){
               reg_allelic <- summary(stats::glm(Disease ~ all_data[, val4] , data = all_data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model
               reg_dominant <- summary(stats::glm(Disease ~  all_data[, val5] , data = all_data, family=stats::binomial(link="logit"), maxit = 100 )) # Logistic regression with dominant model
             }else{
               if(!is.null(covarPath) & is.null(pcaPath)){
                 reg_allelic <- summary(stats::glm(Disease ~ all_data[, val4] + AGE.V0 + SEX + CMV.V0 + NBYTABAC  , data = all_data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model
                 reg_dominant <- summary(stats::glm(Disease ~  all_data[, val5] + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = all_data , family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with dominant model
               }else{
                 if(is.null(covarPath) & !is.null(pcaPath)){
                   reg_allelic <- summary(stats::glm(Disease ~ all_data[, val4] + PC1 + PC2 , data = all_data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model
                   reg_dominant <- summary(stats::glm(Disease ~ all_data[, val5] + PC1 + PC2 , data = all_data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with dominant model


                 }else{
                   reg_allelic <- summary(stats::glm(Disease ~ all_data[, val4] + PC1 + PC2 +  AGE.V0 + SEX + CMV.V0 + NBYTABAC  , data = all_data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model
                   reg_dominant <- summary(stats::glm(Disease ~  all_data[, val5] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = all_data , family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with dominant model
                 }
               }
             }
             #reg_allelic <- summary(stats::glm(Disease ~ all_data[, val4] + PC1 + PC2 + covariables , data = all_data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model
             #reg_dominant <- summary(stats::glm(Disease ~ all_data[, val5] + PC1 + PC2 + covariables, data = all_data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with dominant model

             P_allelic <- reg_allelic$coefficients[2,4]
             OR_allelic <- round(exp(reg_allelic$coefficients[2,1]),2)
             lower_95_CI_allelic <- round(exp(reg_allelic$coefficients[2,1]-1.96*reg_allelic$coefficients[2,2]),2)
             upper_95_CI_allelic <- round(exp(reg_allelic$coefficients[2,1]+1.96*reg_allelic$coefficients[2,2]),2)
             Beta_allelic <- round(reg_allelic$coefficients[2,1],2)
             SE_allelic <- round(reg_allelic$coefficients[2,2],2)
             P_dominant <- reg_dominant$coefficients[2,4]
             OR_dominant <- round(exp(reg_dominant$coefficients[2,1]),2)
             lower_95_CI_dominant <- round(exp(reg_dominant$coefficients[2,1]-1.96*reg_dominant$coefficients[2,2]),2)
             upper_95_CI_dominant <- round(exp(reg_dominant$coefficients[2,1]+1.96*reg_dominant$coefficients[2,2]),2)
             Beta_dominant <- round(reg_dominant$coefficients[2,1],2)
             SE_dominant <- round(reg_dominant$coefficients[2,2],2)
             # add the regression result to the "res" vector
             options(scipen = 999)
             res <- rbind(res, c(hap, freq , freq_control , freq_case, N, P_allelic, OR_allelic, lower_95_CI_allelic, upper_95_CI_allelic, Beta_allelic, SE_allelic, P_dominant, OR_dominant, lower_95_CI_dominant, upper_95_CI_dominant, Beta_dominant, SE_dominant)) # See header of res for elements explanation

           }

           # Save results
           utils::write.table(res,paste(outputFolder , outputFile,"DiseaseAnalysis.csv" ,sep=""), sep=",",row.names=F, quote=F, col.names=F)
           debuggingState(on=FALSE)

           if(file.exists(paste(outputFolder , outputFile,"DiseaseAnalysis.csv" ,sep="")))
           {
             cat("Regression at haplotype levels is Processed successfully \n")
             cat( paste("File", paste(outputFolder , outputFile ,"DiseaseAnalysis.csv" , sep="") ,'is created.  \n', sep = " "))
           }

           ##################################################### PLOT RESULT ############################################################
           information <-   utils::read.table(paste(outputFolder , outputFile,"DiseaseAnalysis.csv" ,sep=""), stringsAsFactors = F, header = T, sep=',')
           seuil <- pvalueCorrection(haploRegResult = paste(outputFolder , outputFile , "-", pheno , ".csv" , sep=""),  correction = "bonferroni")
           if(is.null(seuil)){stop("Seuil calculated with pvalueCorrection() is null")}

           ################### Best Pvalue selection gathered in an unique file
           resultLoop <- as.data.frame(information)
           resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05 ),]
           if(nrow(resultLoop)!=0){
             suppressWarnings({
               options(scipen = 999)
               resultLoop$Phenotype <- pheno
               resultLoop <-resultLoop[,c("Phenotype","HLA_haplotype", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")]
               names(resultLoop) <- c("Phenotype","HLA_haplotype", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")
               write.table(resultLoop, file=paste(outputFolder,"BestPvalue_HaploAnalysis.csv ",sep = "") , append = TRUE, row.names = FALSE ,  col.names = FALSE ,quote=F, sep = ",")
             })
           }

           ################### QQPLOT
           ## Lambda: genomic inflation factor, show the presence of bias
           # allelic Model
           observedR <- sort(information$P_allelic)
           lobsR <- -(log10(observedR))
           expectedR <- c(1:length(observedR))
           lexpR <- -(log10(expectedR / (length(expectedR)+1)))
           lambdaR<-summary(stats::lm(lobsR~lexpR))

           # Dominant Model
           observedD <- sort(information$P_dominant)
           lobsD <- -(log10(observedD))
           expectedD <- c(1:length(observedD))
           lexpD <- -(log10(expectedD / (length(expectedD)+1)))
           lambdaD<- summary(stats::lm(lobsD~lexpD))

           qq1 <- qq(information$P_allelic , main=paste(pheno,"-haplotype (allelic model)",sep = ""), sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
           qq2 <- qq(information$P_dominant , main=paste(pheno,"-haplotype (dominant model)",sep = ""), sub= paste("Lambda = ", round(lambdaD$coefficients[2,1],2),sep = ""))
           capture.output(print(qq1) , file = "NUL")
           capture.output(print(qq2) , file = "NUL")

           ################### MANHATTAN PLOT
           haploSel <- information[information$Frequency >= freq_threshold & (information$P_allelic<=threshold | information$P_dominant<=threshold),]
           if(nrow(haploSel)==0){
             cat(paste("None haplotype to display. You can change your threshold p-value. Set at ",threshold ,sep = ""))
           }else{
             haploSel$order <- c(1:nrow(haploSel))
             options(scipen = 999)
             maxYvalPallelic <- ((max(-log10(haploSel$P_allelic)))+1)
             maxYvalPdom <- ((max(-log10(haploSel$P_dominant)))+1)

             n <- nrow(haploSel)
             if(n<50){ size_point <- 2 }else{ if(n>50 & n<95){ size_point <- 2 }else{ if(n>95){ size_point <- 1}}}

             plot1 <- ggplot( haploSel , aes(x=order, y=  -(log10(P_allelic)) )) + theme(panel.background = element_rect(fill = "#FFE4C4", colour = "#6D9EC1",size = 2, linetype = "solid"))+ geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(haploSel$P_allelic)) > seuil, "red","black")) + scale_x_continuous( limits = c(1 , nrow(haploSel))) + scale_y_continuous( limits = c(0 , maxYvalPallelic)) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")+ ggtitle(paste(pheno, "-haplotype (allelic model)",sep="" ))
             plot2 <- ggplot( haploSel , aes(x=order, y=  -(log10(P_dominant)) )) + theme(panel.background = element_rect(fill = "#FFE4C4", colour = "#6D9EC1",size = 2, linetype = "solid")) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(haploSel$P_dominant)) > seuil, "red","black")) + scale_x_continuous( limits = c(1 ,nrow(haploSel)))  + scale_y_continuous( limits = c(0 , maxYvalPdom)) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")+ ggtitle(paste(pheno, "-haplotype (dominant model)",sep="" ))
             suppressWarnings({
               print(plot1)
               print(plot2)
             })
           }

           grDevices::dev.off()
         }, # end case = disease
         pheno={
           #List of phenotypes
           phenotypes <- stats::na.omit(unique(c(names(phenoFile)))) # Make a list of of phenotypes in the file with colnames
           phenotypes <- phenotypes[phenotypes != "subjectID"]

           # File that will contains best Pvalue for all phenotypes
           X <- c("Phenotype","HLA_haplotype", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")
           write.table(t(X), file=paste(outputFolder,"BestPvalue_HaploAnalysis.csv ",sep = "") , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")

           grDevices::pdf(file = paste(outputFolder ,outputFile,"Haplo_Pheno_analysis_freq",freq_threshold,"_occurency", occurrency ,"_pvalManhattan",threshold , ".pdf" ,sep=""))

           for(pheno in phenotypes)
           {
             print(paste(pheno, "  currently being analysed.", sep = ""))
             ## Merge haplotypes and the phenotype analyzed in this loop
             merge_data_pheno <- merge(all_data,phenoFile[,c("subjectID",pheno)],by = "subjectID")
             colnames(merge_data_pheno)[ncol(merge_data_pheno)] <- "trait"
             phenoGeno <- na.omit(merge_data_pheno)
             Mean_pheno <- mean(phenoGeno$trait)

             ## Determine list of haplotypes to analyze; more than 2 occurences
             df <- as.data.frame(table(stats::na.omit(c(as.character(phenoGeno$A.B.C.DR.DQ.1) , as.character(phenoGeno$A.B.C.DR.DQ.2)))))
             names(df) <- c("haplotypes","Freq")
             haplo <-subset(df , Freq >= occurrency)  ## do a list of haplotypes having at least an "occurrency" times

             if(nrow(haplo)==0){ stop(paste("None haplotype present in your sample appears more than ", occurrency , " times (occurency threshold" , sep=""))}
             res <- c("HLA_haplotype", "Frequency",  "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")

             for(hap in haplo[,1])
             { # Make a loop with all the HLA haplotypes, Convert the haplotypes to 1,2 code (val1, val2), then to genotype (val3), and finally to recessive (val4) or dominant (val5) testable model
               val1 <- paste("allele1_", hap, sep="")
               phenoGeno[,val1] <- ifelse(phenoGeno$A.B.C.DR.DQ.1 == hap, 1,
                                          ifelse(phenoGeno$A.B.C.DR.DQ.1 != hap, 2, NA))

               val2 <- paste("allele2_", hap, sep="")
               phenoGeno[,val2] <- ifelse(phenoGeno$A.B.C.DR.DQ.2 == hap, 1,
                                          ifelse(phenoGeno$A.B.C.DR.DQ.2 != hap, 2, NA))

               val3 <- paste("genotype_", hap, sep="")
               phenoGeno[,val3] <- paste(phenoGeno[,val1], phenoGeno[,val2], sep="")

               val4 <- paste("genoA_", hap, sep="")
               phenoGeno[,val4] <- ifelse(phenoGeno[,val3] == '11', 1,
                                          ifelse(phenoGeno[,val3] == '12', 0,
                                                 ifelse(phenoGeno[,val3] == '21', 0,
                                                        ifelse(phenoGeno[,val3] == '22', -1, NA))))

               val5 <- paste("genoD_", hap, sep="")
               phenoGeno[,val5] <- ifelse(phenoGeno[,val3] == '11', 1,
                                          ifelse(phenoGeno[,val3] == '12', 1,
                                                 ifelse(phenoGeno[,val3] == '21', 1,
                                                        ifelse(phenoGeno[,val3] == '22', 0, NA))))
               ## Frequencies
               freq <- sum(phenoGeno[,val1]==1,phenoGeno[,val2]==1)/(nrow(phenoGeno)*2) # Overall frequency
               n_ind_ctrl <- nrow(subset(phenoGeno, Disease == 0))
               n_ind_cas <-nrow(subset(phenoGeno, Disease == 1))

               if( n_ind_ctrl != 0 & n_ind_cas != 0 )
               { # Frequency are calculated for case and control
                 # N <- 4/((1/nrow(subset(phenoGeno, Disease == 1)))+(1/nrow(subset(phenoGeno, Disease == 0))))
                 N <- nrow(subdata[,!is.na(Disease)])
                 freq_control <- sum(subset(phenoGeno, Disease == 0)[,val1]==1,subset(phenoGeno, Disease == 0)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 0))*2) # Controls frequency
                 freq_case <- sum(subset(phenoGeno, Disease == 1)[,val1]==1,subset(phenoGeno, Disease == 1)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 1))*2) # Cases frequency

               }else{
                 if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl and set to 0
                   freq_control <- 0
                   freq_case <- sum(subset(phenoGeno, Disease == 1)[,val1]==1,subset(phenoGeno, Disease == 1)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 1))*2) # Cases frequency
                   N <- nrow(subset(phenoGeno, Disease == 1))

                 }else{
                   if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for case and set to 0
                     freq_control <- sum(subset(phenoGeno, Disease == 0)[,val1]==1,subset(phenoGeno, Disease == 0)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 0))*2) # Controls frequency
                     freq_case <- 0
                     N <- nrow(subset(phenoGeno, Disease == 0))
                   }else{ # Unknown case control Status
                     N <- 0
                     freq_control <- 0
                     freq_case <- 0
                   }
                 }
               }
               # Adapt the regression correction co-variables
               if(is.null(covarPath) & is.null(pcaPath)){
                 reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, val4] , data = phenoGeno )) # Linear regression with allelic model
                 reg_dominant <- summary(stats::lm(trait ~  phenoGeno[, val5] , data = phenoGeno )) # Linear regression with dominant model
               }else{
                 if(!is.null(covarPath) & is.null(pcaPath)){
                   reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, val4] + AGE.V0 + SEX + CMV.V0 + NBYTABAC  , data = phenoGeno)) # Linear regression with allelic model
                   reg_dominant <- summary(stats::lm(trait ~  phenoGeno[, val5] + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = phenoGeno )) # Linear regression with dominant model
                 }else{
                   if(is.null(covarPath) & !is.null(pcaPath)){
                     reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, val4] + PC1 + PC2, data = phenoGeno)) # Linear regression with allelic model
                     reg_dominant <- summary(stats::lm(trait ~  phenoGeno[, val5] + PC1 + PC2 , data = phenoGeno)) # Linear regression with dominant model
                   }else{
                     reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, val4] + PC1 + PC2 +  AGE.V0 + SEX + CMV.V0 + NBYTABAC  , data = phenoGeno)) # Linear regression with allelic model
                     reg_dominant <- summary(stats::lm(trait ~  phenoGeno[, val5] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = phenoGeno)) # Linear regression with dominant model
                   }
                 }
               }
               # reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, val4] + PC1 + PC2 +  covariables  , data = phenoGeno)) # Logistic regression with allelic model
               # reg_dominant <- summary(stats::lm(trait ~  phenoGeno[, val5] + PC1 + PC2 + covariables , data = phenoGeno )) # Logistic regression with dominant model

               P_allelic <- reg_allelic$coefficients[2,4]
               Beta_allelic <- round(reg_allelic$coefficients[2,1],2)
               SE_allelic <- round(reg_allelic$coefficients[2,2],2)
               Rsquare_allelic <-round(reg_allelic$adj.r.squared, 3)

               P_dominant <- reg_dominant$coefficients[2,4]
               Beta_dominant <- round(reg_dominant$coefficients[2,1],2)
               SE_dominant <- round(reg_dominant$coefficients[2,2],2)
               Rsquare_dominant <-round(reg_dominant$adj.r.squared, 3)

               # add the regression result to the "res" vector
               res <- rbind(res, c(hap, freq , freq_control , freq_case, N, Mean_pheno, P_allelic, Beta_allelic, SE_allelic, Rsquare_allelic, P_dominant, Beta_dominant, SE_dominant, Rsquare_dominant)) # See header of res for elements explanation

             }
             # Save results
             utils::write.table(res, paste(outputFolder , outputFile , "-", pheno , ".csv" ,sep=""), row.names=F, quote=F, col.names=F, sep = ",")

             if(file.exists(paste(outputFolder , outputFile , "-", pheno , ".csv" ,sep="")))
             {
               cat(paste("Regression at haplotype levels with : ", pheno , " expression is processed successfully \n" , sep = ""))
               cat(paste("File", paste(outputFolder , outputFile, "-", pheno , ".csv"  , sep="") ,'is created.  \n', sep = " "))
             }

             ##################################################### PLOT RESULT ############################################################
             information <-   utils::read.table(paste(outputFolder , outputFile , "-", pheno , ".csv" ,sep=""), stringsAsFactors = F, header = T, sep=',')
             seuil <- pvalueCorrection(haploRegResult = paste(outputFolder , outputFile , "-", pheno , ".csv" , sep=""),  correction = "bonferroni")
             if(is.null(seuil)){stop("Seuil calculated with pvalueCorrection() is null")}

             ################### Best Pvalue selection gathered in an unique file
             resultLoop <- as.data.frame(information)
             resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05 ),]
             if(nrow(resultLoop)!=0){
               suppressWarnings({
                 resultLoop$Phenotype <- pheno
                 resultLoop <-resultLoop[,c("Phenotype","HLA_haplotype","Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")]
                 names(resultLoop) <- c("Phenotype","HLA_haplotype", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")
                 write.table(resultLoop, file=paste(outputFolder,"BestPvalue_HaploAnalysis.csv ",sep = "") , append = TRUE, row.names = FALSE ,  col.names = FALSE ,quote=F, sep = ",")
               })
             }

             ################### QQPLOT
             ## Lambda: genomic inflation factor, show the presence of bias
             # allelic Model
             observedR <- sort(information$P_allelic)
             lobsR <- -(log10(observedR))
             expectedR <- c(1:length(observedR))
             lexpR <- -(log10(expectedR / (length(expectedR)+1)))
             lambdaR<-summary(stats::lm(lobsR~lexpR))

             # Dominant Model
             observedD <- sort(information$P_dominant)
             lobsD <- -(log10(observedD))
             expectedD <- c(1:length(observedD))
             lexpD <- -(log10(expectedD / (length(expectedD)+1)))
             lambdaD<- summary(stats::lm(lobsD~lexpD))

             qq1 <- qq(information$P_allelic , main=paste(pheno,"-haplotype (allelic model)",sep = ""), sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
             qq2 <- qq(information$P_dominant , main=paste(pheno,"-haplotype (dominant model)",sep = ""), sub= paste("Lambda = ", round(lambdaD$coefficients[2,1],2),sep = ""))
             capture.output(print(qq1) , file = "NUL")
             capture.output(print(qq2) , file = "NUL")

             ################### MANHATTAN PLOT
             haploSel <- information[information$Frequency >= freq_threshold & (information$P_allelic<=threshold | information$P_dominant<=threshold),]
             if(nrow(haploSel)==0){
               cat(paste("None haplotype to display. You can change your threshold p-value. Set at ",threshold ,sep = ""))
             }else{
               haploSel$order <- c(1:nrow(haploSel))
               options(scipen = 999)
               maxYvalPallelic <- ((max(-log10(haploSel$P_allelic)))+1)
               maxYvalPdom <- ((max(-log10(haploSel$P_dominant)))+1)

               n <- nrow(haploSel)
               if(n<50){ size_point <- 2 }else{ if(n>50 & n<95){ size_point <- 2 }else{ if(n>95){ size_point <- 1}}}

               plot1 <- ggplot( haploSel , aes(x=order, y=  -(log10(P_allelic)) )) + theme(panel.background = element_rect(fill = "#FFE4C4", colour = "#6D9EC1",size = 2, linetype = "solid"))+ geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(haploSel$P_allelic)) > seuil, "red","black")) + scale_x_continuous( limits = c(1 , nrow(haploSel))) + scale_y_continuous( limits = c(0 , maxYvalPallelic)) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")+ ggtitle(paste(pheno, "-haplotype (allelic model)",sep="" ))
               plot2 <- ggplot( haploSel , aes(x=order, y=  -(log10(P_dominant)) )) + theme(panel.background = element_rect(fill = "#FFE4C4", colour = "#6D9EC1",size = 2, linetype = "solid")) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(haploSel$P_dominant)) > seuil, "red","black")) + scale_x_continuous( limits = c(1 ,nrow(haploSel)))  + scale_y_continuous( limits = c(0 , maxYvalPdom)) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")+ ggtitle(paste(pheno, "-haplotype (dominant model)",sep="" ))
               suppressWarnings({
                 print(plot1)
                 print(plot2)
               })
             }

           } # end pheno loop
           grDevices::dev.off()
         }) # end case pheno
}



#######################################################################
###  Amino-acids association analysis with trait(s) (-/+ Covariates)###
#######################################################################

#' Analyze association of HLA Amino-acids with many variate phenotypes (like case/control or quantitative trait)
#' Take as input:
#' (1 required) amino-acids :  subjectID haplotype1 freq1  haplotype2 freq2 postprob followed by many columns of amino-acids("," separator between fields)
#' (2 optional) differents quantitative phenotypes : subjectID traitName1 traitName2 ... traitName(n)
#' (3 required) Disease information : 0= control 1=case (needed for frequencies and/or type=disease analysis)
#' (4 required) the 2 first principals components of the PCA
#' (5 optional) covariable to correct regressions
#'
#' @param aminoAcidFile  HLA Amino-acids data exactly as the one provided by easy-HLA webtool
#' @param phenoPath phenotypes file with different trait names in columns
#' @param statut file conaining subject ID and Disease statut
#' @param pcaPath pca's result: subjectID, PC1 and PC2
#' @param covarPath must contains SubjectID followed by covariables
#' @param outputFolder folder where to save the result
#' @param outputFile result's filename
#' @param type  either "disease" for case/control analysis or "pheno" with quantitative phenotype(s)
#' @return  Linear Regression result and many plots : QQplot and personnalized manhattan plot
#' @import ggplot2
#' @import cowplot
#' @importFrom qqman qq
#' @importFrom utils read.csv write.table capture.output
#' @export


hlaAminoAcidAnalysis  <- function(aminoAcidFile, phenoPath=NULL, statut, pcaPath=NULL, covarPath=NULL, outputFolder ,outputFile, type="pheno")
{
  ## Check parameters
  if(!file.exists(aminoAcidFile)){stop("aminoAcidFile doesn't exist")}
  if(!file.exists(statut)){stop("statut is not found ")}
  if(!is.character(type)){stop("outputFile must be character")}
  if(!dir.exists(outputFolder)){stop("outputFolder directoy doesn't exist. ")}
  if(!is.character(outputFile)){stop("outputFile must be character")}


  ################  DATA IMPORTATION   ################
  amino_acid <- utils::read.csv(aminoAcidFile, stringsAsFactors = F,sep=",", header = TRUE) # Amino acid file from easy-HLA contains : code_patient	A-B-C-DR-DQ-1	freq1	A-B-C-DR-DQ-2	freq2	PostP	resolution FOLLOWED BY Amino acids
  data <- amino_acid[, apply(amino_acid, 2, function(x) length(unique(stats::na.omit(x)))) > 1] # Keep amino acid with information
  colnames(data)[1] <- "subjectID"
  list_aa <- colnames(data)[-c(1:6,(ncol(data)-6):ncol(data))] # list of amino-acids from colnames (first 7 seven colnames are ignored)

  statut <- stats::na.omit(utils::read.table(statut , header = TRUE , sep = "," , dec = "."))
  if(ncol(statut)!= 2){stop("statut must contains 2 colunms : subjectID and Affection statut")}
  names(statut) <- c("subjectID","Disease")
  data <- merge(data, statut, by ="subjectID")

  ################  OPTIONAL PARAMETERS   ################
  if(!is.null(pcaPath)){
    PCs <- stats::na.omit(utils::read.table(pcaPath ,  header = TRUE, sep = "," , dec = "."))
    if(ncol(PCs)!=3){stop("pcaPath should ocntain only 3 columns: subjectID, PC1 and PC2")}
    colnames(PCs)[1]<-"subjectID"
    colnames(PCs)[2]<-"PC1"
    colnames(PCs)[3]<-"PC2"
    data <- merge(data, PCs[,1:3], by ="subjectID")
  }

  if(!is.null(covarPath)){
    covarCorrection <- stats::na.omit(utils::read.table(covarPath ,  header = TRUE, sep = "," , dec = "."))
    colnames(covarCorrection)[1] <- "subjectID"
    covarNames <- stats::na.omit(names(covarCorrection))
    covarNames <- covarNames[covarNames!="subjectID"]
    covariables <- noquote(paste(covarNames, collapse = " + ")) #add as.character ??
    data <- merge(data, covarCorrection, by ="subjectID")
  }


  if(type=="pheno")
  {
    if(is.null(phenoPath)){stop('phenoPath must be provided when type="pheno"')}
    if(!file.exists(phenoPath)){stop("phenoPath file doesn't exist")}
    phenoFile <-  as.data.frame(utils::read.table(phenoPath , stringsAsFactors = F, header = T, sep=','))
    colnames(phenoFile)[1] <- "subjectID"

    #List of phenotypes and covarNames
    phenotypes <- stats::na.omit(unique(c(names(phenoFile)))) # Make a list of of phenotypes in the file with colnames
    phenotypes <- phenotypes[phenotypes != "subjectID"]

  }
  # Thresholds
  freq_threshold <- readline(prompt = "\n \n Enter a frequency threshold of Amino-acids to displays in manhattan plot\n Frequency: " )
  threshold <- readline(prompt = "\n \n In order to display most interesting Amino-acids and in clearly way. \n Define a p value treshold under which one Amino-acids are not displayed in manhattan plot \n P-value:  ")

  switch(type,
         disease={
           res <- c("HLA_Amino_Acid", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic",  "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant","lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")
           X <- c("HLA_Amino_Acid", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")
           write.table(t(X), file=paste(outputFolder,"BestPvalue-AA-DiseaseAssoc.csv ",sep = "") , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")
           grDevices::pdf(file=paste(outputFolder , outputFile , "-AA-DiseaseAssoc_Plot.pdf" ,sep=""))
           for(aa in list_aa) # Make a loop with all amino acids
           {
             ## Frequencies of amino acid
             freq <-  mean(stats::na.omit(data[,aa]))/2 # Overall frequency
             n_ind_ctrl <- nrow(subset(data, Disease == 0))
             n_ind_cas <-nrow(subset(data, Disease == 1))

             if( n_ind_ctrl != 0 & n_ind_cas != 0 )
             { # Frequency are calculated for case and control
               N <- nrow(data[,!is.na(Disease)])
               freq_control <- mean(stats::na.omit(subset(data, Disease == 0)[,aa]))/2 # Controls frequency
               freq_case <- mean(stats::na.omit(subset(data, Disease == 1)[,aa]))/2 # Cases frequency

             }else{
               if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl and set to 0
                 freq_control <- 0
                 freq_case <- mean(stats::na.omit(subset(data, Disease == 1)[,aa]))/2 # Cases frequency
                 N <- nrow(subset(data, Disease == 1))

               }else{
                 if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for case and set to 0
                   freq_control <-  mean(stats::na.omit(subset(data, Disease == 0)[,aa]))/2 # Controls frequency
                   freq_case <- 0
                   N <- nrow(subset(data, Disease == 0))
                 }else{ # Unknown case control Status
                   stop("For a Case/Control analysis, the Disease colunm is required. Please note 1=case and 0=control")
                 }
               }
             }


             if(is.null(covarPath) & is.null(pcaPath)){
               # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
               data[,aa][data[,aa] == 0] <- -1
               data[,aa][data[,aa] == 1] <- 0
               data[,aa][data[,aa] == 2] <- 1
               reg_allelic <- summary(stats::glm(Disease ~ data[, aa] , data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

               # Convert (-1,0,1) to (0,1) to run the regression in dominant model
               data[,aa][data[,aa] == 0] <- 1
               data[,aa][data[,aa] == -1] <- 0
               reg_dominant <- summary(stats::glm(Disease ~ data[, aa] , data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

             }else{
               if(!is.null(covarPath) & is.null(pcaPath)){
                 # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                 data[,aa][data[,aa] == 0] <- -1
                 data[,aa][data[,aa] == 1] <- 0
                 data[,aa][data[,aa] == 2] <- 1
                 reg_allelic <- summary(stats::glm(Disease ~ data[, aa] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                 # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                 data[,aa][data[,aa] == 0] <- 1
                 data[,aa][data[,aa] == -1] <- 0
                 reg_dominant <- summary(stats::glm(Disease ~ data[, aa] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model
               }else{
                 if(is.null(covarPath) & !is.null(pcaPath)){
                   # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                   data[,aa][data[,aa] == 0] <- -1
                   data[,aa][data[,aa] == 1] <- 0
                   data[,aa][data[,aa] == 2] <- 1
                   reg_allelic <- summary(stats::glm(Disease ~ data[, aa] + PC1 + PC2, data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                   # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                   data[,aa][data[,aa] == 0] <- 1
                   data[,aa][data[,aa] == -1] <- 0
                   reg_dominant <- summary(stats::glm(Disease ~ data[, aa] + PC1 + PC2, data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                 }else{
                   # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                   data[,aa][data[,aa] == 0] <- -1
                   data[,aa][data[,aa] == 1] <- 0
                   data[,aa][data[,aa] == 2] <- 1
                   reg_allelic <- summary(stats::glm(Disease ~ data[, aa] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                   # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                   data[,aa][data[,aa] == 0] <- 1
                   data[,aa][data[,aa] == -1] <- 0
                   reg_dominant <- summary(stats::glm(Disease ~ data[, aa] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = data, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                 }}}
             # add the regression result to the "res" vector
             res <- rbind(res, c(aa, freq, freq_control, freq_case,N ,reg_allelic$coefficients[2,4], exp(reg_allelic$coefficients[2,1]), exp(reg_allelic$coefficients[2,1]-1.96*reg_allelic$coefficients[2,2]), exp(reg_allelic$coefficients[2,1]+1.96*reg_allelic$coefficients[2,2]), reg_allelic$coefficients[2,1], reg_allelic$coefficients[2,2], reg_dominant$coefficients[2,4], exp(reg_dominant$coefficients[2,1]), exp(reg_dominant$coefficients[2,1]-1.96*reg_dominant$coefficients[2,2]), exp(reg_dominant$coefficients[2,1]+1.96*reg_dominant$coefficients[2,2]), reg_dominant$coefficients[2,1], reg_dominant$coefficients[2,2])) # See header of res for elements explanation

           }
           grDevices::dev.off()
           # Save results
           utils::write.table(res, paste(outputFolder , outputFile,"-AA-DiseaseAssoc-Analysis" ,".csv", sep=""), row.names=F, quote=F, col.names=F, sep = ",")
           if(file.exists( paste(outputFolder , outputFile,"-AA-DiseaseAssoc-Analysis" ,".csv", sep=""))){
             cat("Regression at Amino-Acid levels is Processed successfully \n")
             cat(paste("File", paste(outputFolder , outputFile , sep="") ,'is created. \n', sep = " "))
           }
           ## Import Files
           aaReg_result <- utils::read.table(paste(outputFolder , outputFile,"-AA-DiseaseAssoc-Analysis" ,".csv", sep="") , header= TRUE , sep = "," , dec = "." , na.strings = NA)
           seuil <- pvalueCorrection(AARegResult = paste(outputFolder , outputFile,"-AA-DiseaseAssoc-Analysis" ,".csv", sep=""), correction = "bonferroni")
           if(is.null(seuil)){stop("Seuil calculated with pvalueCorrection() is null")}

           # Ne garder que la position et AA
           aaReg_result$AA <- gsub("^[A-Z]{2}_[A-Z]_","",aaReg_result$HLA_Amino_Acid )
           aaReg_result$AA <- gsub("^[A-Z]{2}_[A-Z]{1,3}\\d_","",aaReg_result$AA)

           aaReg_result$Locus <- gsub("^[A-Z][A-Z]_", "", aaReg_result$HLA_Amino_Acid) # keep gene's name
           aaReg_result$Locus <- gsub("_\\d{1,}_[A-Z]", "",  aaReg_result$Locus)
           aaReg_result$Locus <- gsub("_[[:punct:]]\\d{1,}_[A-Z]", "",  aaReg_result$Locus)
           aaReg_result$Locus <- gsub("_\\d{1,}_[[:punct:]]", "",  aaReg_result$Locus)

           ################### Best Pvalue selection gathered in an unique file
           resultLoop <- as.data.frame(aaReg_result)
           resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold  & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05) ,]
           if(nrow(resultLoop)!=0){
             suppressWarnings({
               resultLoop$Phenotype <- pheno
               resultLoop <-resultLoop[,c("Phenotype","Locus","AA", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")]
               write.table(resultLoop, file=paste(outputFolder,"BestPvalue-AA-DiseaseAssoc.csv ",sep = "") , append = TRUE, row.names = FALSE ,  col.names = FALSE ,quote=F, sep = ",")
             })

           }
           ################### QQPLOT
           ## Lambda: genomic inflation factor, show the presence of bias
           # allelic Model
           observedR <- na.omit(sort(aaReg_result$P_allelic))
           lobsR <- -(log10(observedR))
           expectedR <- c(1:length(observedR))
           lexpR <- -(log10(expectedR / (length(expectedR)+1)))
           lambdaR<-summary(lm(lobsR~lexpR))

           # Dominant Model
           observedD <- sort(aaReg_result$P_dominant)
           lobsD <- -(log10(observedD))
           expectedD <- c(1:length(observedD))
           lexpD <- -(log10(expectedD / (length(expectedD)+1)))
           lambdaD<- summary(lm(lobsD~lexpD))

           ## QQ PLOT
           qq1 <- qq(aaReg_result$P_allelic ,main="Allelic model Disease-Amino-acids assoc" ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
           qq2 <- qq(aaReg_result$P_dominant , main="Dominant model Disease-Amino-acids assoc" ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
           capture.output(print(qq1),file='NUL')
           capture.output(print(qq2),file='NUL')

           ## Graphic for visualizing P value of AA
           aaSel <- aaReg_result[aaReg_result$Frequency >= freq_threshold & (aaReg_result$P_allelic<=threshold | aaReg_result$P_dominant<=threshold),]
           aSel <- aaSel[, c("AA","Locus","P_allelic","P_dominant")]
           if(nrow(aaSel)==0){ message(paste("None Amino-Acid to display. You can change your frequency and threshold p-value.\n Set at :\n Frequency : ", freq_threshold, "\n Pvalue threshold :",threshold,sep = ""))}
           plot5 <- ggplot( aaSel, aes(x=AA, y=  -(log10(P_allelic)) )) + geom_point() + facet_grid(.~Locus , scales = "free" , space = "free_x") + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red") + ggtitle(" Amino-Acid Analysis in allelic Model")
           plot6 <- ggplot( aaSel , aes(x=AA, y=  -(log10(P_dominant)) )) + geom_point() + facet_grid(.~Locus , scales = "free" , space = "free_x") + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")+ ggtitle(" Amino-Acid Analysis in Dominant Model")

           suppressWarnings({
             print(plot5)
             print(plot6)
           })

           grDevices::dev.off()
         }, # end case = disease
         pheno ={

           # File that will contains best Pvalue for all phenotypes
           X <- c("Phenotype","Locus","AA", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")
           write.table(t(X), file=paste(outputFolder,"BestPvalue-AA-PhenoAssoc.csv",sep = "") , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")

           grDevices::pdf(file=paste(outputFolder , outputFile , "-AA-PhenotypeAssoc_regressionPlot.pdf" ,sep=""),10,10)
           for(pheno in phenotypes)
           {
             res <- c("HLA_Amino_Acid", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")
             print(paste(pheno, " is currently being analysed.", sep = ""))
             ## Merge AA and the phenotype analyzed in this loop
             merge_data_pheno <- merge(data,phenoFile[,c("subjectID",pheno)],by = "subjectID")
             colnames(merge_data_pheno)[ncol(merge_data_pheno)] <- "trait"
             phenoGeno <- na.omit(merge_data_pheno)
             Mean_pheno <- mean(phenoGeno$trait)

             for(aa in list_aa) # Make a loop with all amino acids
             {
               # Caculate frequency of amino acid
               freq <- mean(stats::na.omit(phenoGeno[,aa]))/2

               n_ind_ctrl <- nrow(subset(phenoGeno, Disease == 0))
               n_ind_cas <-nrow(subset(phenoGeno, Disease == 1))

               if( n_ind_ctrl != 0 & n_ind_cas != 0 )
               { # Frequency are calculated for case and control
                 N <- 4/((1/nrow(subset(phenoGeno, Disease == 1)))+(1/nrow(subset(phenoGeno, Disease == 0))))
                 freq_control <- mean(stats::na.omit(subset(phenoGeno, Disease == 0)[,aa]))/2
                 freq_case <- mean(stats::na.omit(subset(phenoGeno, Disease == 1)[,aa]))/2

               }else{
                 if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl
                   freq_control <- 0
                   freq_case <-mean(stats::na.omit(subset(phenoGeno, Disease == 1)[,aa]))/2
                   N <- nrow(subset(phenoGeno, Disease == 1))

                 }else{
                   if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for case
                     # N <- 4/(1/nrow(subset(subdata, Disease == 0)))
                     freq_control <- mean(stats::na.omit(subset(phenoGeno, Disease == 0)[,aa]))/2
                     freq_case <- 0
                     N <- nrow(subset(phenoGeno, Disease == 0))

                   }else{ # Unknown case control Status. This info is not essential for the analysis
                     N <- 0
                     freq_control <- 0
                     freq_case <- 0
                   }}}

               # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
               # phenoGeno[,aa][phenoGeno[,aa] == 0] <- -1
               # phenoGeno[,aa][phenoGeno[,aa] == 1] <- 0
               # phenoGeno[,aa][phenoGeno[,aa] == 2] <- 1

               # genoA <- paste("genoA_",aa,sep="")
               # phenoGeno[,genoA] <- ifelse(phenoGeno[,aa][phenoGeno[,aa] == 0],-1,
               # ifelse(phenoGeno[,aa][phenoGeno[,aa] == 1],0,
               # ifelse(phenoGeno[,aa][phenoGeno[,aa] == 2],1, NA)))
               # Convert (-1,0,1) to (0,1) to run the regression in dominant model
               # phenoGeno[,aa][phenoGeno[,aa] == 0] <- 1
               # phenoGeno[,aa][phenoGeno[,aa] == -1] <- 0
               # genoD <- paste("genoD_",aa,sep="")
               # phenoGeno[, genoD] <- ifelse(phenoGeno[,aa][phenoGeno[,aa] == 0],1,
               #                                ifelse(phenoGeno[,aa][phenoGeno[,aa] == -1],0,NA))
               # print(phenoGeno[,genoA])
               # print(phenoGeno[,genoD])
               if(is.null(covarPath) & is.null(pcaPath)){

                 reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, aa] , data = phenoGeno)) # Linear regression with allelic model
                 reg_dominant <- summary(stats::lm(trait ~ phenoGeno[, aa] , data = phenoGeno)) # Linear regression with allelic model

               }else{
                 if(!is.null(covarPath) & is.null(pcaPath)){
                   # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                   phenoGeno[,aa][phenoGeno[,aa] == 0] <- -1
                   phenoGeno[,aa][phenoGeno[,aa] == 1] <- 0
                   phenoGeno[,aa][phenoGeno[,aa] == 2] <- 1
                   reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, aa] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno )) # Linear regression with allelic model

                   # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                   phenoGeno[,aa][phenoGeno[,aa] == 0] <- 1
                   phenoGeno[,aa][phenoGeno[,aa] == -1] <- 0
                   phenoGeno[, aa] <-  as.factor(phenoGeno[, aa])
                   reg_dominant <- summary(stats::lm(trait ~ data[, aa] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno )) # Linear regression with allelic model
                 }else{
                   if(is.null(covarPath) & !is.null(pcaPath)){
                     # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                     phenoGeno[,aa][phenoGeno[,aa] == 0] <- -1
                     phenoGeno[,aa][phenoGeno[,aa] == 1] <- 0
                     phenoGeno[,aa][phenoGeno[,aa] == 2] <- 1
                     reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, aa] + PC1 + PC2, data = phenoGeno)) # Linear regression with allelic model

                     # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                     phenoGeno[,aa][phenoGeno[,aa] == 0] <- 1
                     phenoGeno[,aa][phenoGeno[,aa] == -1] <- 0
                     reg_dominant <- summary(stats::lm(trait ~ phenoGeno[, aa] + PC1 + PC2, data = phenoGeno)) # Linear regression with allelic model

                   }else{
                     # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                     phenoGeno[,aa][phenoGeno[,aa] == 0] <- -1
                     phenoGeno[,aa][phenoGeno[,aa] == 1] <- 0
                     phenoGeno[,aa][phenoGeno[,aa] == 2] <- 1
                     reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, aa] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno)) # Linear regression with allelic model

                     # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                     phenoGeno[,aa][phenoGeno[,aa] == 0] <- 1
                     phenoGeno[,aa][phenoGeno[,aa] == -1] <- 0
                     reg_dominant <- summary(stats::lm(trait ~ phenoGeno[, aa] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno )) # Linear regression with allelic model

                   }}}

               P_allelic <-  reg_allelic$coefficients[2,4]
               Beta_allelic <- round(reg_allelic$coefficients[2,1],3)
               SE_allelic <- round(reg_allelic$coefficients[2,2],3)
               Rsquare_allelic <-round(reg_allelic$adj.r.squared, 3)

               P_dominant <- reg_dominant$coefficients[2,4]
               Beta_dominant <- round(reg_dominant$coefficients[2,1],3)
               SE_dominant <- round(reg_dominant$coefficients[2,2],3)
               Rsquare_dominant <-round(reg_dominant$adj.r.squared, 3)

               # add the regression result to the "res" vector
               res <- rbind(res, c(aa,  freq, freq_control, freq_case, N, Mean_pheno, P_allelic, Beta_allelic, SE_allelic, Rsquare_allelic, P_dominant, Beta_dominant, SE_dominant, Rsquare_dominant)) # See header of res for elements explanation

             }

             # Save results
             utils::write.table(res, paste(outputFolder , outputFile,"-AA-", pheno,"-Analysis" ,".csv", sep=""), row.names=F, quote=F, col.names=F, sep = ",")
             if(file.exists( paste(outputFolder , outputFile,"-AApheno-Analysis" ,".csv", sep=""))){
               cat(paste("Regression at amino-acid levels with : ", pheno , " expression is processed successfully \n" , sep = ""))
               cat(paste("File", paste(outputFolder , outputFile , sep="") ,'is created. \n', sep = " "))
             }

             ## Import Files
             aaReg_result <- utils::read.table( paste(outputFolder , outputFile,"-AA-", pheno,"-Analysis" ,".csv", sep="") , header= TRUE , sep = "," , dec = "." , na.strings = NA)
             seuil <- pvalueCorrection(AARegResult =  paste(outputFolder , outputFile,"-AA-", pheno,"-Analysis" ,".csv", sep=""), correction = "bonferroni")
             if(is.null(seuil)){stop("Bonferroni threshold is null ! Line 1333 ")}
             # Ne garder que la position et AA
             aaReg_result$AA <- gsub("^[A-Z]{2}_[A-Z]_","",aaReg_result$HLA_Amino_Acid )
             aaReg_result$AA <- gsub("^[A-Z]{2}_[A-Z]{1,3}\\d_","",aaReg_result$AA)
             aaReg_result$Locus <- gsub("^[A-Z][A-Z]_", "", aaReg_result$HLA_Amino_Acid) # keep gene's name
             aaReg_result$Locus <- gsub("_\\d{1,}_[A-Z]", "",  aaReg_result$Locus)
             aaReg_result$Locus <- gsub("_[[:punct:]]\\d{1,}_[A-Z]", "",  aaReg_result$Locus)
             aaReg_result$Locus <- gsub("_\\d{1,}_[[:punct:]]", "",  aaReg_result$Locus)
             aaReg_result$Locus <- as.factor(aaReg_result$Locus)

             ################### Best Pvalue selection gathered in an unique file
             resultLoop <- as.data.frame(aaReg_result)
             resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold  & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05) ,]
             if(nrow(resultLoop)!=0){
               suppressWarnings({
                 resultLoop$Phenotype <- pheno
                 resultLoop <-resultLoop[,c("Phenotype","Locus","AA","Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")]
                 write.table(resultLoop, file=paste(outputFolder,"BestPvalue-AA-PhenoAssoc.csv",sep = ""), append = TRUE , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")
               })
             }
             ################### QQPLOT
             ## Lambda: genomic inflation factor, show the presence of bias
             # allelic Model
             observedR <- na.omit(sort(aaReg_result$P_allelic))
             lobsR <- -(log10(observedR))
             expectedR <- c(1:length(observedR))
             lexpR <- -(log10(expectedR / (length(expectedR)+1)))
             lambdaR<-summary(lm(lobsR~lexpR))

             # Dominant Model
             observedD <- sort(aaReg_result$P_dominant)
             lobsD <- -(log10(observedD))
             expectedD <- c(1:length(observedD))
             lexpD <- -(log10(expectedD / (length(expectedD)+1)))
             lambdaD<- summary(lm(lobsD~lexpD))

             ## QQ PLOT
             qq1 <- qq(aaReg_result$P_allelic ,main=paste(" AA vs ",pheno," (allelic model)",sep = "") ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
             qq2 <- qq(aaReg_result$P_dominant , main=paste(" AA vs ",pheno," (dominant model)",sep = "") ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
             capture.output(print(qq1),file='NUL')
             capture.output(print(qq2),file='NUL')

             ## Graphic for visualizing P value of AA
             aaSel <- aaReg_result[aaReg_result$Frequency >= freq_threshold & (aaReg_result$P_allelic<=threshold | aaReg_result$P_dominant<=threshold),]
             aSel <- aaSel[, c("AA","Locus","P_allelic","P_dominant")]
             if(nrow(aaSel)==0){ message(paste("None Amino-Acid to display. You can change your frequency and threshold p-value.\n Set at :\n Frequency : ", freq_threshold, "\n Pvalue threshold :",threshold,sep = ""))}

             ######### FONCTIONNE MAIS ILLISIBLE AXE X
             ## Adapt x axis legend with the number of alleles to plot
             n <- nrow(aaSel)
             if(n<50){
               size_txt_x_axis <- 10
               size_point <- 3
             }else{
               if(n>50 & n<95){
                 size_txt_x_axis <- 8
                 size_point <- 3
               }else{
                 if(n>95 & n<130){
                   size_txt_x_axis <- 7
                   size_point <- 2
                 }else{
                   if(n>130){
                     size_txt_x_axis <- 6
                     size_point <- 2
                   }}}}

             panelNames <- list( # change ggplot panel's names
               "A" = "HLA-A",
               "B" = "HLA-B ",
               "C" = "HLA-C",
               "DRB1" = "HLA-DRB1 ",
               "DQB1" = "HLA-DQB1"
             )
             variable_labeller <- function(variable,value){
               return(panelNames[value])
             }

             ## allelic model
             aaA <-  aaSel[aaSel$Locus=="A",]
             plot1 <- ggplot( aaA, aes(x=AA, y=  -(log10(P_allelic)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaA$P_allelic)) > seuil, "red","black")) + facet_grid( .~ Locus  ) +geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")+ xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaB <-  aaSel[aaSel$Locus=="B",]
             plot2 <- ggplot( aaB, aes(x=AA, y=  -(log10(P_allelic)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaB$P_allelic)) > seuil, "red","black"))  + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red") + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis-1, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaC <-  aaSel[aaSel$Locus=="C",]
             plot3 <- ggplot( aaC, aes(x=AA, y=  -(log10(P_allelic)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaC$P_allelic)) > seuil, "red","black")) + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red") + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis-1, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaDR <-  aaSel[aaSel$Locus=="DRB1",]
             plot4 <- ggplot( aaDR, aes(x=AA, y=  -(log10(P_allelic)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaDR$P_allelic)) > seuil, "red","black")) + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red") + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaDQ <-  aaSel[aaSel$Locus=="DQB1",]
             plot5 <- ggplot( aaDQ, aes(x=AA, y=  -(log10(P_allelic)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaDQ$P_allelic)) > seuil, "red","black")) + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")  + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             plotclass1 <- ggdraw() +
               draw_label(paste(" AA vs ",pheno," (allelic model)",sep = ""), fontface='bold',x = 0.5, y = .97, hjust = 0.5, vjust = 0.5,color="#ec6649",size=12) +
               draw_plot(plot1, 0, .5, 1, .45) +
               draw_plot(plot2, 0, 0, .5, .5) +
               draw_plot(plot3, .5, 0, .5, .5)
             plotclass2 <- ggdraw() +
               draw_label(paste(" AA vs ",pheno," (allelic model)",sep = ""), fontface='bold',x = 0.5, y =  .97, hjust = 0.5, vjust = 0.5,color="#ec6649",size=12) +
               draw_plot(plot4, 0, .5, 1, .45) +
               draw_plot(plot5, 0, 0, 1, .5)
             suppressWarnings({
               print(plotclass1)
               print(plotclass2)
             })

             ## dominant model
             aaA <-  aaSel[aaSel$Locus=="A",]
             plot1 <- ggplot( aaA, aes(x=AA, y=  -(log10(P_dominant)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaA$P_dominant)) > seuil, "red","black")) + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")+ xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaB <-  aaSel[aaSel$Locus=="B",]
             plot2 <- ggplot( aaB, aes(x=AA, y=  -(log10(P_dominant)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaB$P_dominant)) > seuil, "red","black"))  + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red") + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis-1, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaC <-  aaSel[aaSel$Locus=="C",]
             plot3 <- ggplot( aaC, aes(x=AA, y=  -(log10(P_dominant)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaC$P_dominant)) > seuil, "red","black")) + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red") + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis-1, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaDR <-  aaSel[aaSel$Locus=="DRB1",]
             plot4 <- ggplot( aaDR, aes(x=AA, y=  -(log10(P_dominant)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaDR$P_dominant)) > seuil, "red","black")) + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red") + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             aaDQ <-  aaSel[aaSel$Locus=="DQB1",]
             plot5 <- ggplot( aaDQ, aes(x=AA, y=  -(log10(P_dominant)) )) + geom_point(na.rm = TRUE, size=size_point, color=ifelse((-log10(aaDQ$P_dominant)) > seuil, "red","black")) + facet_grid( .~ Locus  ) + geom_hline(yintercept= -log10(0.05) , color = "blue") + geom_hline(yintercept= seuil , color = "red")  + xlab("position and amino-acid")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis, angle=90), axis.title.y = element_text(color="#2b79be", size=11),strip.text = element_text(face="bold",size = 10), strip.background = element_rect(fill="lightblue"))
             # pour changer en rouge le texte de l'axe x comme les points, mais ne foncionne pas
             #theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(size = size_txt_x_axis, angle=90, color=ifelse((-log10(aaDQ$P_dominant)) > seuil, "red","black"))
             plotclass1 <- ggdraw() +
               draw_label(paste(" AA vs ",pheno," (dominant model)",sep = ""), fontface='bold',x = 0.5, y = .97, hjust = 0.5, vjust = 0.5,color="#ec6649",size=12) +
               draw_plot(plot1, 0, .5, 1, .45) +
               draw_plot(plot2, 0, 0, .5, .5) +
               draw_plot(plot3, .5, 0, .5, .5)
             plotclass2 <- ggdraw() +
               draw_label(paste(" AA vs ",pheno," (dominant model)",sep = ""), fontface='bold',x = 0.5, y =  .97, hjust = 0.5, vjust = 0.5,color="#ec6649",size=12) +
               draw_plot(plot4, 0, .5, 1, .45) +
               draw_plot(plot5, 0, 0, 1, .5)
             suppressWarnings({
               print(plotclass1)
               print(plotclass2)
             })

           }})# end case = pheno
  grDevices::dev.off()
}




############################################################################
####  KIR ligands association analysis with trait(s) (-/+ Covariates)   ####
############################################################################

#' Analyze association of HLA KIR ligands with many variate phenotypes (like case/control or quantitative trait)
#' Take as input:
#' (1 required) KIR ligands  :  subjectID haplotype1 freq1  haplotype2 freq2 postprob ("," separator between fields)
#' (2 optional) differents quantitative phenotypes : subjectID traitName1 traitName2 ... traitName(n)
#' (3 required) Disease information : 0= control 1=case (needed for frequencies and/or type=disease analysis)
#' (4 required) the 2 first principals components of the PCA
#' (5 optional) covariable to correct regressions
#'
#' @param datapath HLA File from Easy-HLA
#' @param phenoPath phenotypes file with different trait names in columns
#' @param statut file conaining subject ID and Disease statut
#' @param pcaPath pca result for PC1 and PC2
#' @param covarPath must contains SubjectID followed by covariables
#' @param outputFolder Folder's pathways where to save the result file
#' @param outputFile Name of the result file
#' @param type etiher "disease" for analysis in case/control or "pheno" for analysis with an quantitative phenotype
#' @import ggplot2
#' @importFrom qqman qq
#' @importFrom utils read.csv write.csv capture.output
#' @return Regression result in a file with columns: Amino_Acid	Frequency	Frequency_control	Frequency_case	N	P_allelic	OR_allelic	lower_95%_CI_allelic	upper_95%_CI_allelic	Beta_allelic	SE_allelic	P_dominant	OR_dominant	lower_95%_CI_dominant	upper_95%_CI_dominant	Beta_dominant	SE_dominant
#' @export


hlaKIRAnalysis <- function(datapath, phenoPath=NULL, statut, pcaPath=NULL, covarPath=NULL, outputFolder ,outputFile, type="pheno")
{

  ## Check parameters
  if(is.null(datapath)){stop("datapath must be provided")}
  if(!file.exists(datapath)){stop("datapath doesn't exist")}
  if(is.null(statut)){stop("statut must be provided of both analysis type")}

  ## Import data
  # Import hlaKIR file from easy-HLA code_patient	A-B-C-DR-DQ-1	freq1	A-B-C-DR-DQ-2	freq2	PostP	resolution 3DL2 Bw4 Bw4I Bw4T Bw6 C1 C2
  kirLigand <- utils::read.csv(datapath, stringsAsFactors = F,sep=",", header = TRUE) # Import kir ligand acid data
  data <- kirLigand[, apply(kirLigand, 2, function(x) length(unique(stats::na.omit(x)))) > 1] # Keep kir ligand colums with information
  if(ncol(data) < 7){stop("datapath must contains more than 7 colunms like the file provided by easy-HLA :code_patient,	A-B-C-DR-DQ-1,	freq1,	A-B-C-DR-DQ-2,	freq2,	PostP,	resolution, followed by kir ligand. For example (3DL2, Bw4, Bw4I, Bw4T, Bw6, C1, C2)")}

  colnames(data)[1] <- "subjectID"
  list <- colnames(data)
  list_ligand <- list[7:length(list)]

  ## Merge
  # Optional parameters
  if(!is.null(pcaPath)){
    if(!file.exists(pcaPath)){stop("pcaPath doesn't exist")}
    PCs <- stats::na.omit(utils::read.table(pcaPath ,  header = TRUE, sep = "," , dec = "."))
    if(ncol(PCs)!=3){stop("pcaPath should ocntain only 3 columns: subjectID, PC1 and PC2")}
    colnames(PCs)[1]<-"subjectID"
    colnames(PCs)[2]<-"PC1"
    colnames(PCs)[3]<-"PC2"
    data <- merge(data, PCs[,1:3], by ="subjectID")
  }

  if(!is.null(covarPath)){
    if(!file.exists(covarPath)){stop("covarPath doesn't exist")}
    covarCorrection <- stats::na.omit(utils::read.table(covarPath ,  header = TRUE, sep = "," , dec = "."))
    colnames(covarCorrection)[1] <- "subjectID"
    covarNames <- stats::na.omit(names(covarCorrection))
    covarNames <- covarNames[covarNames!="subjectID"]
    covariables <- noquote(paste(covarNames, collapse = " + ")) #add as.character ??
    data <- merge(data, covarCorrection, by ="subjectID")
  }

  if(type=="disease")
  {
    statut <- stats::na.omit(utils::read.table(statut , header = TRUE , sep = "," , dec = "."))
    if(ncol(statut)!= 2){stop("statut must contains 2 colunms : subjectID and Affection statut")}
    names(statut) <- c("subjectID","Disease")
    data <- merge(data, statut, by ="subjectID")

  }else{
    if(type=="pheno")
    {
      if(is.null(phenoPath)){stop('phenoPath must be provided when type="pheno"')}
      if(!file.exists(phenoPath)){stop("phenoPath file doesn't exist")}
      phenoFile <-  as.data.frame(utils::read.table(phenoPath , stringsAsFactors = F, header = T, sep=','))
      colnames(phenoFile)[1] <- "subjectID"
      # List of phenotypes and covarNames
      phenotypes <- stats::na.omit(unique(c(names(phenoFile)))) # Make a list of of phenotypes in the file with colnames
      phenotypes <- phenotypes[phenotypes != "subjectID"]

      statut <- stats::na.omit(utils::read.table(statut , header = TRUE , sep = "," , dec = "."))

      if(ncol(statut)!=2){stop("statut mut contains 2 colunms : subjectID and Affection statut")}
      names(statut) <- c("subjectID","Disease")
      data <- merge(data, statut, by ="subjectID")

    }
  }
  switch(type,
         disease={
           res <- c("HLA_KIR", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic",  "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant","lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")
           freq_threshold <- readline(prompt = "\n \n Enter a frequency threshold of KIR-ligand to displays in manhattan plot\n Frequency: " )
           threshold <- readline(prompt = "\n \n In order to display most interesting KIR-ligand and in clearly way. \n Define a p value treshold under which one haplotypes are not displayed in manhattan plot \n P-value:  ")
           X <- c("HLA_KIR", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")
           write.table(t(X), file=paste(outputFolder,"BestPvalue-kir-DisAssoc.csv",sep = "") , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")

           pdf(file=paste(outputFolder , outputFile , "-KIR-DiseaseAssocPlot.pdf" ,sep=""),10,10)

           ## Merge kir and the phenotype analyzed in this loop
           phenoGeno <- na.omit(data)

           for(kir in list_ligand) # Make a loop with all kir ligands
           {
             # Caculate frequency of the kir ligand
             freq <- mean(stats::na.omit(phenoGeno[, kir]))/2

             n_ind_ctrl <- nrow(subset(phenoGeno, Disease == 0))
             n_ind_cas <-nrow(subset(phenoGeno, Disease == 1))

             if( n_ind_ctrl != 0 & n_ind_cas != 0 )
             { # Frequency are calculated for case and control
               N <- 4/((1/nrow(subset(phenoGeno, Disease == 1)))+(1/nrow(subset(phenoGeno, Disease == 0))))
               freq_control <- mean(stats::na.omit(subset(phenoGeno, Disease == 0)[,kir]))/2
               freq_case <- mean(stats::na.omit(subset(phenoGeno, Disease == 1)[,kir]))/2

             }else{
               if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl
                 freq_control <- 0
                 freq_case <-mean(stats::na.omit(subset(phenoGeno, Disease == 1)[,kir]))/2
                 N <- nrow(subset(phenoGeno, Disease == 1))

               }else{
                 if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for case
                   # N <- 4/(1/nrow(subset(subdata, Disease == 0)))
                   freq_control <- mean(stats::na.omit(subset(phenoGeno, Disease == 0)[,kir]))/2
                   freq_case <- 0
                   N <- nrow(subset(phenoGeno, Disease == 0))

                 }else{ # Unknown case control Status
                   N <- 0
                   freq_control <- 0
                   freq_case <- 0
                 }}}

             if(is.null(covarPath) & is.null(pcaPath)){
               # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
               phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
               phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
               phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
               reg_allelic <- summary(stats::glm(Disease ~ phenoGeno[, kir] , data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

               # Convert (-1,0,1) to (0,1) to run the regression in dominant model
               phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
               phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
               reg_dominant <- summary(stats::glm(Disease ~ phenoGeno[, kir] , data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

             }else{
               if(!is.null(covarPath) & is.null(pcaPath)){
                 # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                 phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
                 phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
                 phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
                 reg_allelic <- summary(stats::glm(Disease ~ phenoGeno[, kir] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100 )) # Logistic regression with allelic model

                 # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                 phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
                 phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
                 phenoGeno[, kir] <-  as.factor(phenoGeno[, kir])
                 reg_dominant <- summary(stats::glm(Disease ~ data[, kir] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno , family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model
               }else{
                 if(is.null(covarPath) & !is.null(pcaPath)){
                   # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                   phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
                   phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
                   phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
                   reg_allelic <- summary(stats::glm(Disease ~ phenoGeno[, kir] + PC1 + PC2, data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                   # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                   phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
                   phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
                   reg_dominant <- summary(stats::glm(Disease ~ phenoGeno[, kir] + PC1 + PC2, data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                 }else{
                   # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                   phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
                   phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
                   phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
                   reg_allelic <- summary(stats::glm(Disease ~ phenoGeno[, kir] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)) # Logistic regression with allelic model

                   # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                   phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
                   phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
                   reg_dominant <- summary(stats::glm(Disease ~ phenoGeno[, kir] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno )) # Logistic regression with allelic model

                 }}}

             P_allelic <-  reg_allelic$coefficients[2,4]
             OR_allelic <-  round(exp(reg_allelic$coefficients[2,1]),3)
             lower_95_CI_allelic <- round(exp(reg_allelic$coefficients[2,1]-1.96*reg_allelic$coefficients[2,2]),3)
             upper_95_CI_allelic <- round(exp(reg_allelic$coefficients[2,1]+1.96*reg_allelic$coefficients[2,2]),3)
             Beta_allelic <- round(reg_allelic$coefficients[2,1],3)
             SE_allelic <- round(reg_allelic$coefficients[2,2],3)
             P_dominant <- reg_dominant$coefficients[2,4]
             OR_dominant <- round(exp(reg_dominant$coefficients[2,1]),3)
             lower_95_CI_dominant <- round(exp(reg_dominant$coefficients[2,1]-1.96*reg_dominant$coefficients[2,2]),3)
             upper_95_CI_dominant <- round(exp(reg_dominant$coefficients[2,1]+1.96*reg_dominant$coefficients[2,2]),3)
             Beta_dominant <- round(reg_dominant$coefficients[2,1],3)
             SE_dominant <- round(reg_dominant$coefficients[2,2],3)
             options(scipen = 999)
             # add the regression result to the "res" vector
             res <- rbind(res, c(kir, freq, freq_control, freq_case, N ,P_allelic,OR_allelic,lower_95_CI_allelic,upper_95_CI_allelic,Beta_allelic,SE_allelic,P_dominant,OR_dominant,lower_95_CI_dominant,upper_95_CI_dominant,Beta_dominant,SE_dominant))
           }

           # Save results
           utils::write.table(res, paste(outputFolder , outputFile,"-kir-Disease-Analysis" ,".csv", sep=""), row.names=F, quote=F, col.names=F, sep = ",")
           if(file.exists( paste(outputFolder , outputFile,"-kir-Disease-Analysis" ,".csv", sep=""))){
             cat("Regression of KIR ligands is Processed successfully \n")
             cat(paste("File", paste(outputFolder , outputFile , sep="") ,'is created. \n', sep = " "))
           }

           ## Import Files
           kirReg_result <- utils::read.table( paste(outputFolder , outputFile,"-kir-Disease-Analysis" ,".csv", sep="") , header= TRUE , sep = "," , dec = "." , na.strings = NA)
           seuil <- pvalueCorrection(KIR =  paste(outputFolder , outputFile,"-kir-Disease-Analysis" ,".csv", sep=""), correction = "bonferroni")
           if(is.null(seuil)){stop("Bonferroni threshold is null ! ")}

           ################### Best Pvalue selection gathered in an unique file
           resultLoop <- as.data.frame(kirReg_result)
           resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold  & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05) ,]
           # if(resultLoop!=factor(0)){ # if ==factor(0) error is displayed "argument is of length zero"
           if(!is.null(resultLoop) & nrow(resultLoop)!=0){
             suppressWarnings({
               resultLoop$Phenotype <- pheno
               resultLoop <-resultLoop[,c("HLA_KIR", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_0.95_CI_dominant", "upper_0.95_CI_dominant", "Beta_dominant", "SE_dominant")]
               write.table(resultLoop, file=paste(outputFolder,"BestPvalue-kir-DisAssoc.csv",sep = ""), append = TRUE , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")
             })
           }
           # }
           ################### QQPLOT
           ## Lambda: genomic inflation factor, show the presence of bias
           # allelic Model
           observedR <- na.omit(sort(kirReg_result$P_allelic))
           lobsR <- -(log10(observedR))
           expectedR <- c(1:length(observedR))
           lexpR <- -(log10(expectedR / (length(expectedR)+1)))
           lambdaR<-summary(lm(lobsR~lexpR))

           # Dominant Model
           observedD <- sort(kirReg_result$P_dominant)
           lobsD <- -(log10(observedD))
           expectedD <- c(1:length(observedD))
           lexpD <- -(log10(expectedD / (length(expectedD)+1)))
           lambdaD<- summary(lm(lobsD~lexpD))

           ## QQ PLOT
           qq1 <- qq(kirReg_result$P_allelic ,main=paste(" kir vs ",pheno," (allelic model)",sep = "") ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
           qq2 <- qq(kirReg_result$P_dominant , main=paste(" kir vs ",pheno," (dominant model)",sep = "") ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
           capture.output(print(qq1),file='NUL')
           capture.output(print(qq2),file='NUL')

           ## Graphic for visualizing P value of kir
           kirSel <- kirReg_result[kirReg_result$Frequency >= freq_threshold & (kirReg_result$P_allelic<=threshold | kirReg_result$P_dominant<=threshold),]
           kirSel <- kirSel[, c("HLA_KIR","P_allelic","P_dominant")]
           if(nrow(kirSel)==0){ message(paste("None KIR ligand to display for ", pheno , ". You can change your frequency and threshold p-value.\n Set at :\n Frequency : ", freq_threshold, "\n Pvalue threshold :",threshold,sep = ""))}

           ######### FONCTIONNE MAIS ILLISIBLE AXE X
           ## Adapt x axis legend with the number of alleles to plot


           ## allelic model
           plot1 <- ggplot( kirSel, aes(x=HLA_KIR, y=  -(log10(P_allelic)))) + geom_point(na.rm = TRUE, color=ifelse((-log10(kirSel$P_allelic)) >= seuil, "red","black")) + geom_hline(yintercept= -log10(0.05) , color = "blue") + ggtitle(paste(" kir vs ",pheno," (allelic model)",sep = ""))+  geom_hline(yintercept= seuil , color = "red")+ xlab("KIR ligand")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(angle=0), axis.title.y = element_text(color="#2b79be", size=11))
           ## dominant model
           plot2 <- ggplot( kirSel, aes(x=HLA_KIR, y=  -(log10(P_dominant)))) + geom_point(na.rm = TRUE, color=ifelse((-log10(kirSel$P_dominant)) >= seuil, "red","black")) + geom_hline(yintercept= -log10(0.05) , color = "blue") + ggtitle(paste(" kir vs ",pheno," (dominant model)",sep = ""))+  geom_hline(yintercept= seuil , color = "red")+ xlab("KIR ligand")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(angle=0), axis.title.y = element_text(color="#2b79be", size=11))
           suppressWarnings({
             print(plot1)
             print(plot2)
           })


           dev.off()

         },
         pheno ={
           freq_threshold <- readline(prompt = "\n \n Enter a frequency threshold of KIR-ligand to displays in manhattan plot\n Frequency: " )
           threshold <- readline(prompt = "\n \n In order to display most interesting KIR-ligand and in clearly way. \n Define a p value treshold under which one haplotypes are not displayed in manhattan plot \n P-value:  ")
           X <- c("Phenotype","HLA_KIR","Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")
           write.table(t(X), file=paste(outputFolder,"BestPvalue-kir-PhenoAssoc.csv",sep = "") , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")

           pdf(file=paste(outputFolder , outputFile , "-KIR-PhenotypeAssocPlot.pdf" ,sep=""),10,10)
           for(pheno in phenotypes)
           {
             # pheno <- "MFI_CD16_in_CD16hi_of_NKnew.panel4"
             res <- c("HLA_KIR", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")

             ## Merge kir and the phenotype analyzed in this loop
             merge_data_pheno <- merge(data,phenoFile[,c("subjectID",pheno)],by = "subjectID")
             colnames(merge_data_pheno)[ncol(merge_data_pheno)] <- "trait"
             phenoGeno <- na.omit(merge_data_pheno)

             for(kir in list_ligand) # Make a loop with all amino acids
             {
               # Caculate frequency of amino acid
               freq <- mean(stats::na.omit(phenoGeno[, kir]))/2

               n_ind_ctrl <- nrow(subset(phenoGeno, Disease == 0))
               n_ind_cas <-nrow(subset(phenoGeno, Disease == 1))

               if( n_ind_ctrl != 0 & n_ind_cas != 0 )
               { # Frequency are calculated for case and control
                 N <- 4/((1/nrow(subset(phenoGeno, Disease == 1)))+(1/nrow(subset(phenoGeno, Disease == 0))))
                 freq_control <- mean(stats::na.omit(subset(phenoGeno, Disease == 0)[,kir]))/2
                 freq_case <- mean(stats::na.omit(subset(phenoGeno, Disease == 1)[,kir]))/2

               }else{
                 if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl
                   freq_control <- 0
                   freq_case <-mean(stats::na.omit(subset(phenoGeno, Disease == 1)[,kir]))/2
                   N <- nrow(subset(phenoGeno, Disease == 1))

                 }else{
                   if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for case
                     # N <- 4/(1/nrow(subset(subdata, Disease == 0)))
                     freq_control <- mean(stats::na.omit(subset(phenoGeno, Disease == 0)[,kir]))/2
                     freq_case <- 0
                     N <- nrow(subset(phenoGeno, Disease == 0))

                   }else{ # Unknown case control Status
                     N <- 0
                     freq_control <- 0
                     freq_case <- 0
                   }}}

               if(is.null(covarPath) & is.null(pcaPath)){
                 # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                 phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
                 phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
                 phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
                 reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, kir] , data = phenoGeno)) # Linear regression with allelic model

                 # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                 phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
                 phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
                 reg_dominant <- summary(stats::lm(trait ~ phenoGeno[, kir] , data = phenoGeno)) # Linear regression with allelic model

               }else{
                 if(!is.null(covarPath) & is.null(pcaPath)){
                   # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                   phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
                   phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
                   phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
                   reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, kir] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno )) # Linear regression with allelic model

                   # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                   phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
                   phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
                   phenoGeno[, kir] <-  as.factor(phenoGeno[, kir])
                   reg_dominant <- summary(stats::lm(trait ~ data[, kir] + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno )) # Linear regression with allelic model
                 }else{
                   if(is.null(covarPath) & !is.null(pcaPath)){
                     # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                     phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
                     phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
                     phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
                     reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, kir] + PC1 + PC2, data = phenoGeno)) # Linear regression with allelic model

                     # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                     phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
                     phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
                     reg_dominant <- summary(stats::lm(trait ~ phenoGeno[, kir] + PC1 + PC2, data = phenoGeno)) # Linear regression with allelic model

                   }else{
                     # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
                     phenoGeno[,kir][phenoGeno[,kir] == 0] <- -1
                     phenoGeno[,kir][phenoGeno[,kir] == 1] <- 0
                     phenoGeno[,kir][phenoGeno[,kir] == 2] <- 1
                     reg_allelic <- summary(stats::lm(trait ~ phenoGeno[, kir] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno)) # Linear regression with allelic model

                     # Convert (-1,0,1) to (0,1) to run the regression in dominant model
                     phenoGeno[,kir][phenoGeno[,kir] == 0] <- 1
                     phenoGeno[,kir][phenoGeno[,kir] == -1] <- 0
                     reg_dominant <- summary(stats::lm(trait ~ phenoGeno[, kir] + PC1 + PC2 + AGE.V0 + SEX + CMV.V0 + NBYTABAC, data = phenoGeno )) # Linear regression with allelic model

                   }}}

               P_allelic <-  reg_allelic$coefficients[2,4]
               Beta_allelic <- round(reg_allelic$coefficients[2,1],3)
               SE_allelic <- round(reg_allelic$coefficients[2,2],3)
               Rsquare_allelic <-round(reg_allelic$adj.r.squared, 3)

               P_dominant <- reg_dominant$coefficients[2,4]
               Beta_dominant <- round(reg_dominant$coefficients[2,1],3)
               SE_dominant <- round(reg_dominant$coefficients[2,2],3)
               Rsquare_dominant <-round(reg_dominant$adj.r.squared, 3)

               options(scipen = 999)

               res <- rbind(res, c(kir, freq, freq_control, freq_case, N, Mean_pheno, P_allelic, Beta_allelic, SE_allelic, Rsquare_allelic, P_dominant, Beta_dominant, SE_dominant, Rsquare_dominant)) # See header of res for elements explanation

             }

             # Save results
             utils::write.table(res, paste(outputFolder , outputFile,"-kir-", pheno,"-Analysis" ,".csv", sep=""), row.names=F, quote=F, col.names=F, sep = ",")
             if(file.exists( paste(outputFolder , outputFile,"-kirpheno-Analysis" ,".csv", sep=""))){
               cat("Regression on KIR ligands is Processed successfully \n")
               cat(paste("File", paste(outputFolder , outputFile , sep="") ,'is created. \n', sep = " "))
             }

             ## Import Files
             kirReg_result <- utils::read.table( paste(outputFolder , outputFile,"-kir-", pheno,"-Analysis" ,".csv", sep="") , header= TRUE , sep = "," , dec = "." , na.strings = NA)
             seuil <- pvalueCorrection(kirRegResult =  paste(outputFolder , outputFile,"-kir-", pheno,"-Analysis" ,".csv", sep=""), correction = "bonferroni")
             if(is.null(seuil)){stop("Bonferroni threshold is null ! ")}

             ################### Best Pvalue selection gathered in an unique file
             resultLoop <- as.data.frame(kirReg_result)
             resultLoop <- resultLoop[resultLoop$Frequency >= freq_threshold  & (resultLoop$P_allelic <= 0.05 | resultLoop$P_dominant <= 0.05) ,]
             # if(resultLoop!=factor(0)){ # if ==factor(0) error is displayed "argument is of length zero"
             if(!is.null(resultLoop) & nrow(resultLoop)!=0){
               suppressWarnings({
                 resultLoop$Phenotype <- pheno
                 resultLoop <-resultLoop[,c("Phenotype","HLA_KIR", "Frequency", "Frequency_control", "Frequency_case", "N","Mean_pheno", "P_allelic", "Beta_allelic", "SE_allelic", "R2_allelic","P_dominant", "Beta_dominant", "SE_dominant","R2_dominant")]
                 write.table(resultLoop, file=paste(outputFolder,"BestPvalue-kir-PhenoAssoc.csv",sep = ""), append = TRUE , row.names = FALSE ,   col.names = FALSE  ,quote=F, sep = ",")
               })
             }
             # }
             ################### QQPLOT
             ## Lambda: genomic inflation factor, show the presence of bias
             # allelic Model
             observedR <- na.omit(sort(kirReg_result$P_allelic))
             lobsR <- -(log10(observedR))
             expectedR <- c(1:length(observedR))
             lexpR <- -(log10(expectedR / (length(expectedR)+1)))
             lambdaR<-summary(lm(lobsR~lexpR))

             # Dominant Model
             observedD <- sort(kirReg_result$P_dominant)
             lobsD <- -(log10(observedD))
             expectedD <- c(1:length(observedD))
             lexpD <- -(log10(expectedD / (length(expectedD)+1)))
             lambdaD<- summary(lm(lobsD~lexpD))

             ## QQ PLOT
             qq1 <- qq(kirReg_result$P_allelic ,main=paste(" kir vs ",pheno," (allelic model)",sep = "") ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
             qq2 <- qq(kirReg_result$P_dominant , main=paste(" kir vs ",pheno," (dominant model)",sep = "") ,  sub= paste("Lambda = ", round(lambdaR$coefficients[2,1],2),sep = ""))
             capture.output(print(qq1),file='NUL')
             capture.output(print(qq2),file='NUL')

             ## Graphic for visualizing P value of kir
             kirSel <- kirReg_result[kirReg_result$Frequency >= freq_threshold & (kirReg_result$P_allelic<=threshold | kirReg_result$P_dominant<=threshold),]
             kirSel <- kirSel[, c("HLA_KIR","P_allelic","P_dominant")]
             if(nrow(kirSel)==0){ message(paste("None KIR ligand to display for ", pheno , ". You can change your frequency and threshold p-value.\n Set at :\n Frequency : ", freq_threshold, "\n Pvalue threshold :",threshold,sep = ""))}

             ######### FONCTIONNE MAIS ILLISIBLE AXE X
             ## Adapt x axis legend with the number of alleles to plot


             ## allelic model
             plot1 <- ggplot( kirSel, aes(x=HLA_KIR, y=  -(log10(P_allelic)))) + geom_point(na.rm = TRUE, color=ifelse((-log10(kirSel$P_allelic)) >= seuil, "red","black")) + geom_hline(yintercept= -log10(0.05) , color = "blue") + ggtitle(paste(" kir vs ",pheno," (allelic model)",sep = ""))+  geom_hline(yintercept= seuil , color = "red")+ xlab("KIR ligand")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(angle=0), axis.title.y = element_text(color="#2b79be", size=11))
             ## dominant model
             plot2 <- ggplot( kirSel, aes(x=HLA_KIR, y=  -(log10(P_dominant)))) + geom_point(na.rm = TRUE, color=ifelse((-log10(kirSel$P_dominant)) >= seuil, "red","black")) + geom_hline(yintercept= -log10(0.05) , color = "blue") + ggtitle(paste(" kir vs ",pheno," (dominant model)",sep = ""))+  geom_hline(yintercept= seuil , color = "red")+ xlab("KIR ligand")  + theme(axis.title.x = element_text(color="#2b79be", size=11), axis.text.x = element_text(angle=0), axis.title.y = element_text(color="#2b79be", size=11))

             suppressWarnings({
               print(plot1)
               print(plot2)
             })

           }
           dev.off()
         })


}


#######################################################
###  HLA-C expression analysis with many phenotypes ###
#######################################################


#' HLA-C expression analys many variate phenotypes
#'
#' @param dataPath file with subject.ID and HLA-c expression "quantitative"
#' @param phenoPath  file with subject.ID and as many columns as quantitative phenotypes
#' @param statut file conaining subject ID and Disease statut
#' @param pcaPath pca result for PC1 and PC2
#' @param covarPath must contains SubjectID followed by covariables
#' @param outputFolder folder where to save the result
#' @param outputFile Name for the result csv file with columns: subjectID	Disease	PC1	PC2	A	A	B	B	C	C	DRB1	DRB1	DQB1	DQB1
#' @param type etiher "disease" for analysis in case/control or "pheno" for analysis with an quantitative phenotype
#' @return Linear Regression result with the phenotypes or case/control and corrected either with PCs or with covariates
#' @import ggplot2
#' @importFrom utils read.table write.table capture.output
#' @importFrom qqman qq
#' @importFrom stats binomial glm na.omit p.adjust lm
#' @importFrom grDevices dev.off  pdf
#' @export
#'

hlaCexpressionAnalysis <- function( dataPath , phenoPath=NULL  , statut=NULL, pcaPath=NULL,  covarPath=NULL , outputFolder , outputFile , type="pheno")
{
  # Obligatory
  data <-  as.data.frame(utils::read.table(dataPath , stringsAsFactors = F, header = T, sep=','))
  if(ncol(data)!=2){stop("datapath should ocntain only 2 columns: subjectID and HLAC expression")}
  names(data) <- c("subjectID" , "HLAexprC")

  # Option parameters
  if(!is.null(pcaPath)){
    PCs <- stats::na.omit(utils::read.table(pcaPath ,  header = TRUE, sep = "," , dec = "."))
    if(ncol(PCs)!=3){stop("pcaPath should ocntain only 3 columns: subjectID, PC1 and PC2")}
    colnames(PCs)[1]<-"subjectID"
    colnames(PCs)[2]<-"PC1"
    colnames(PCs)[3]<-"PC2"
    data <- merge(data, PCs[,1:3], by ="subjectID")
  }

  if(!is.null(covarPath)){
    covarCorrection <- stats::na.omit(utils::read.table(covarPath ,  header = TRUE, sep = "," , dec = "."))
    colnames(covarCorrection)[1] <- "subjectID"
    covarNames <- stats::na.omit(names(covarCorrection))
    covarNames <- covarNames[covarNames!="subjectID"]
    covariables <- noquote(paste(covarNames, collapse = " + ")) #add as.character ??
    data <- merge(data, covarCorrection, by ="subjectID")
  }

  if(type=="disease")
  {
    if(is.null(statut)){stop("statut must be provided")}
    if(!file.exists(statut)){stop("statut doesn't exist")}
    statut <- stats::na.omit(utils::read.table(statut , header = TRUE , sep = "," , dec = "."))
    if(ncol(statut)!=2){stop("statut mut contains 2 colunms : subjectID and Affection statut")}
    names(statut) <- c("subjectID","Disease")
    data <- merge(data, statut, by ="subjectID")

  }else{
    if(type=="pheno")
    {
      if(is.null(phenoPath)){stop('phenoPath must be provided when type="pheno"')}
      if(!file.exists(phenoPath)){stop("phenoPath file doesn't exist")}
      phenoFile <-  as.data.frame(utils::read.table(phenoPath , stringsAsFactors = F, header = T, sep=','))
      colnames(phenoFile)[1] <- "subjectID"
      #List of phenotypes and covarNames
      phenotypes <- stats::na.omit(unique(c(names(phenoFile)))) # Make a list of of phenotypes in the file with colnames
      phenotypes <- phenotypes[phenotypes != "subjectID"]

    }
  }
  ggplotRegression <- function (fit) {
    requireNamespace(ggplot2)
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      labs(title = paste(pheno," vs HLA-C Expression"),
           subtitle=paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                          "Intercept =",signif(fit$coef[[1]],5 ),
                          " Slope =",signif(fit$coef[[2]], 5),
                          " P =",signif(summary(fit)$coef[2,4], 5))) +
      xlab("HLA-C expression ") + ylab(paste(pheno," expression", sep=""))+
      theme(
        axis.title.x = element_text(color="#2b79be", size=11),
        axis.title.y = element_text(color="#2b79be", size=11),
        plot.title = element_text(hjust = 0.5, size = 12, color = "#ec6649", face = "bold"),    # Center title position and size
        plot.subtitle = element_text(hjust = 0.5,size=10, color = "#07985b",face = "italic"))            # Center subtitle
  }


  switch(type,
         disease={
           res <- c("Mean_HLACexpr","N","N_ctrl","N_case", "P", "Beta", "SE","R2")
           pdf(file=paste(outputFolder , outputFile , "-HLA-KIR_Disease_analysis_regressionPlot.pdf" ,sep=""))

           ## Frequencies
           freq <- sum(phenoGeno[,val1]==1,phenoGeno[,val2]==1)/(nrow(phenoGeno)*2) # Overall frequency
           n_ind_ctrl <- nrow(subset(phenoGeno, Disease == 0))
           n_ind_cas <-nrow(subset(phenoGeno, Disease == 1))

           if( n_ind_ctrl != 0 & n_ind_cas != 0 )
           { # Frequency are calculated for case and control
             # N <- 4/((1/nrow(subset(all_data, Disease == 1)))+(1/nrow(subset(all_data, Disease == 0))))
             N <- nrow(phenoGeno[,!is.na(Disease)])
             freq_control <- sum(subset(phenoGeno, Disease == 0)[,val1]==1,subset(phenoGeno, Disease == 0)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 0))*2) # Controls frequency
             freq_case <- sum(subset(phenoGeno, Disease == 1)[,val1]==1,subset(phenoGeno, Disease == 1)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 1))*2) # Cases frequency

           }else{
             if( n_ind_ctrl == 0 & n_ind_cas != 0){ # won't calculate frequency for ctrl and set to 0
               freq_control <- 0
               freq_case <- sum(subset(phenoGeno, Disease == 1)[,val1]==1,subset(phenoGeno, Disease == 1)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 1))*2) # Cases frequency
               N <- nrow(subset(phenoGeno, Disease == 1))

             }else{
               if( n_ind_ctrl != 0 & n_ind_cas == 0){ # won't calculate frequency for case and set to 0
                 freq_control <- sum(subset(phenoGeno, Disease == 0)[,val1]==1,subset(phenoGeno, Disease == 0)[,val2]==1)/(nrow(subset(phenoGeno, Disease == 0))*2) # Controls frequency
                 freq_case <- 0
                 N <- nrow(subset(phenoGeno, Disease == 0))
               }else{ # Unknown case control Status
                 N <- 0
                 freq_control <- 0
                 freq_case <- 0
               }
             }
           }
           if(is.null(covarPath) & is.null(pcaPath)){
             reg <- summary(stats::glm(Disease~HLAexprC, data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)) # explain the phenotype with HLAC expression
           }else{
             if(!is.null(covarPath) & is.null(pcaPath)){
               reg <- summary(stats::glm(Disease~HLAexprC + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)) # explain the phenotype with HLAC expression
             }else{
               if(is.null(covarPath) & !is.null(pcaPath)){
                 reg <- stats::glm(Disease~HLAexprC + PC1 + PC2 , data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100)
               }else{
                 reg <- summary(stats::glm(Disease~HLAexprC + AGE.V0 + SEX + CMV.V0 + NBYTABAC + PC1 + PC2 , data = phenoGeno, family=stats::binomial(link="logit"), maxit = 100))
               }}}

           # print(reg)
           P <- reg$coefficients[2,4]
           OR <- round(exp(reg$coefficients[2,1]),2)
           lower_95_CI <- round(exp(reg$coefficients[2,1]-1.96*reg$coefficients[2,2]),2)
           upper_95_CI <- round(exp(reg$coefficients[2,1]+1.96*reg$coefficients[2,2]),2)
           Beta <- round(reg$coefficients[2,1],2)
           SE <- round(reg$coefficients[2,2],2)
           MeanHLACexpr <- mean(phenoGeno$HLAexprC)
           Rsquare <-signif(reg$adj.r.squared, 5)
           # add the regression result to the "res" vector
           res <- rbind(res, c(MeanHLACexpr, N,Nctrl, Ncase, P, OR,Rsquare, lower_95_CI, upper_95_CI, Beta, SE)) # See header of res for elements explanation

           print(ggplotRegression(reg))
           dev.off()
           utils::write.table(res, paste(outputFolder , outputFile , "-HLAC-expr-Analysis" , ".csv" ,sep=""), row.names=F, quote=F, col.names=F, sep = ",")
           if(file.exists(paste(outputFolder , outputFile , "-HLAC-expr-Analysis" , ".csv" ,sep=""))){
             cat(paste("File", paste(outputFolder , outputFile, "-", pheno , ".csv"  , sep="") ,'is created.  \n', sep = " "))
           }
         },
         pheno ={
           res <- c("Phenotype","Mean_phenotype","Mean_HLACexpr","N", "P", "OR","R^2", "lower_0.95_CI", "upper_0.95_CI", "Beta", "SE")
           pdf(file=paste(outputFolder , outputFile , "-regressionPlot.pdf" ,sep=""))

           for(pheno in phenotypes)
           {
             # Merge of pheno columns and alleles
             phenoGeno <- merge(data,phenoFile[,c("subjectID",pheno)],by = "subjectID")
             colnames(phenoGeno)[ncol(phenoGeno)] <- "trait"
             phenoGeno <- na.omit(phenoGeno)

             N <- nrow(phenoGeno)
             if(is.null(covarPath) & is.null(pcaPath)){
               reg <- summary(stats::lm(trait~HLAexprC, data = phenoGeno)) # explain the phenotype with HLAC expression

             }else{
               if(!is.null(covarPath) & is.null(pcaPath)){
                 reg <- summary(stats::lm(trait~HLAexprC + AGE.V0 + SEX + CMV.V0 + NBYTABAC , data = phenoGeno)) # explain the phenotype with HLAC expression
               }else{
                 if(is.null(covarPath) & !is.null(pcaPath)){
                   reg <- stats::lm(trait~HLAexprC + PC1 + PC2 , data = phenoGeno)
                 }else{reg <- summary(stats::lm(trait~HLAexprC + AGE.V0 + SEX + CMV.V0 + NBYTABAC + PC1 + PC2 , data = phenoGeno))}}}

             P <- summary(reg)$coefficients[2,4]
             lower_95_CI <- round(exp(summary(reg)$coefficients[2,1]-1.96*summary(reg)$coefficients[2,2]),2)
             upper_95_CI <- round(exp(summary(reg)$coefficients[2,1]+1.96*summary(reg)$coefficients[2,2]),2)
             Beta <- round(summary(reg)$coefficients[2,1],2)
             SE <- round(summary(reg)$coefficients[2,2],2)
             Mean_pheno <- mean(phenoGeno$trait)
             MeanHLACexpr <- mean(phenoGeno$HLAexprC)
             Rsquare <-signif(summary(reg)$adj.r.squared, 5)
             # add the regression result to the "res" vector
             res <- rbind(res, c(pheno,Mean_pheno,MeanHLACexpr, N, P,Rsquare, lower_95_CI, upper_95_CI, Beta, SE)) # See header of res for elements explanation

             print(ggplotRegression(reg))

           }
           dev.off()
           utils::write.table(res, paste(outputFolder , outputFile , "-HLAC-expr-Analysis" , ".csv" ,sep=""), row.names=F, quote=F, col.names=F, sep = ",")
           if(file.exists(paste(outputFolder , outputFile , "-HLAC-expr-Analysis" , ".csv" ,sep=""))){
             cat("Regression of HLA-C expression is Processed successfully \n")
             cat(paste("File", paste(outputFolder , outputFile , "-HLAC-expr-Analysis" , ".csv"  , sep="") ,'is created.  \n', sep = " "))
           }})

}


#' Function that correct P-value for HLA alleles, haplotypes or Amino Acids, kir ligand and HLA-C expression. Propose 2 differents methods: bonferroni with allele's regression result or Permutation by redoing regression over 1000 times or different number of loop.
#' @param allele (permutation) Imputation's result file with 14 colunms : subjectID, Disease, PC1, PC2, prob, A, A, B, B, C, C, DRB1, DRB1, DQB1, DQB1 ("," separator between fields)
#' @param haplo (permutation) Haplotypes file with at least these columns : subjectID,	A-B-C-DR-DQ-1,	freq1,	A-B-C-DR-DQ-2,	freq2,	PostP ("," separator between fields)
#' @param AA  (permutation) Amino-Acid File from easy-hla
#' @param KIR (permutation) KIR ligand's file from easy-hla
#' @param HLAC (permutation) HLAC expression as provided by easy-hla
#' @param alleleRegResult (bonferroni) result file of Regression on allele : function
#' @param haploRegResult (bonferroni) result file of Regression on haplotypes : function
#' @param AARegResult (bonferroni) result file of Regression on AA : function
#' @param kirRegResult (bonferroni) result file of Regression on kir ligands : function
#' @param hlaCRegResult (bonferroni) result file of Regression on kir ligands : function
#' @param pcaFile FIle with sample id and two principal components of PCA
#' @param statut File with sample id and disease statut (coded 0 1)
#' @param correction Method to correct P value: permutation or bonferroni (default).
#' @return Return a list with the new critical P value for each level of analysis: seuilAllele, seuilHaplo and seuilAA
#' @importFrom utils read.table write.table capture.output
#' @export
#'

pvalueCorrection <- function(allele = NULL, haplo = NULL, AA = NULL, KIR = NULL, HLAC= NULL
                             ,alleleRegResult = NULL, haploRegResult = NULL, AARegResult = NULL, kirRegResult=NULL, hlaCRegResult=NULL
                             ,pcaFile = NULL , statut = NULL , correction = "bonferroni"){
  # 1er choix la correction
  # 2eme choix le niveau d'analyse
  # Calcul d'un seul seuil puis retourner par la fonction
  # #' @param numberLoop Number of loop for regression when performing permutation correction of P value.
  # , numberLoop=100
  if( correction == "bonferroni" ){
    if(!is.null(alleleRegResult)){
      data <- utils::read.table(alleleRegResult , header= TRUE , sep = "," , dec = "." , na.strings = NA)
      seuil <- -log10((0.05/(nrow(data)-1))) # nombre d'allèles uniques exsitants
    }
    if(!is.null(haploRegResult)){
      data <- utils::read.table(haploRegResult , header= TRUE , sep = "," , dec = "." , na.strings = NA)
      seuil <- -log10((0.05/(nrow(data)-1)))

    }
    if(!is.null(AARegResult)){
      data <- utils::read.table(AARegResult , header= TRUE , sep = "," , dec = "." , na.strings = NA)
      seuil <- -log10((0.05/(nrow(data)-1)))
    }
    if(!is.null(kirRegResult)){
      data <- utils::read.table(kirRegResult , header= TRUE , sep = "," , dec = "." , na.strings = NA)
      seuil <- -log10((0.05/(nrow(data)-1)))
    }
    if(!is.null(hlaCRegResult)){
      data <- utils::read.table(hlaCRegResult , header= TRUE , sep = "," , dec = "." , na.strings = NA)
      seuil <- -log10((0.05/(nrow(data)-1)))
    }

  }else{

    if( correction == "permutation"){

      if(!is.null(allele)){
        data <- utils::read.table(allele , header= TRUE , sep = "," , dec = "." , na.strings = NA)
        HLA_genes <- c("A", "B", "C", "DRB1", "DQB1")
        res <- c("HLA_allele", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_95%_CI_allelic", "upper_95%_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_95%_CI_dominant", "upper_95%_CI_dominant", "Beta_dominant", "SE_dominant")
        p <- c()
        best_Psample <- c()
        for( i in 1:1000) # set to 1000 by default
        {
          print(i)
          pheno <- sample(data$Disease, replace = TRUE )
          data$Disease <- pheno

          for(gene in HLA_genes) # Make a loop with all the HLA genes
          {
            subdata <- data[,c("subjectID","Disease","PC1","PC2",gene,paste(gene,".1",sep = ""))]
            names(subdata) <- c("subjectID","Disease","PC1","PC2", "allele1" , "allele2" )
            HLA <- na.omit(as.character(unique(c(unique(subdata$allele1), unique(subdata$allele2))))) # Make a list of all existing HLA alleles in each file

            for(allele in HLA)
            {
              val1 <- paste("allele1_", allele, sep="")
              subdata[,val1] <- ifelse(subdata$allele1 == allele, 1,
                                       ifelse(subdata$allele1 != allele, 2, NA))
              val2 <- paste("allele2_", allele, sep="")
              subdata[,val2] <- ifelse(subdata$allele2 == allele, 1,
                                       ifelse(subdata$allele2 != allele, 2, NA))
              val3 <- paste("genotype_", allele, sep="")
              subdata[,val3] <- paste(subdata[,val1], subdata[,val2], sep="")
              val4 <- paste("genoA_", allele, sep="")
              subdata[,val4] <- ifelse(subdata[,val3] == '11', 1,
                                       ifelse(subdata[,val3] == '12', 0,
                                              ifelse(subdata[,val3] == '21', 0,
                                                     ifelse(subdata[,val3] == '22', -1, NA))))
              val5 <- paste("genoD_", allele, sep="")
              subdata[,val5] <- ifelse(subdata[,val3] == '11', 1,
                                       ifelse(subdata[,val3] == '12', 1,
                                              ifelse(subdata[,val3] == '21', 1,
                                                     ifelse(subdata[,val3] == '22', 0, NA))))
              reg_allelic <- summary(glm(Disease ~ subdata[, val4] + PC1 + PC2, data = subdata, family=binomial(link="logit"))) # Logistic regression with allelic model
              reg_dominant <- summary(glm(Disease ~ subdata[, val5] + PC1 + PC2, data = subdata, family=binomial(link="logit"))) # Logistic regression with dominant model

              freq <- sum(subdata[,val1]==1,subdata[,val2]==1)/(nrow(subdata)*2) # Overall frequency
              freq_control <- sum(subset(subdata, Disease == 0)[,val1]==1,subset(subdata, Disease == 0)[,val2]==1)/(nrow(subset(subdata, Disease == 0))*2) # Controls frequency
              freq_case <- sum(subset(subdata, Disease == 1)[,val1]==1,subset(subdata, Disease == 1)[,val2]==1)/(nrow(subset(subdata, Disease == 1))*2) # Cases frequency with homozygote for the allele (11)
              res <- rbind(res, c(allele , freq , freq_control , freq_case, 4/((1/nrow(subset(data, Disease == 1)))+(1/nrow(subset(data, Disease == 0)))), reg_allelic$coefficients[2,4], exp(reg_allelic$coefficients[2,1]), exp(reg_allelic$coefficients[2,1]-1.96*reg_allelic$coefficients[2,2]), exp(reg_allelic$coefficients[2,1]+1.96*reg_allelic$coefficients[2,2]), reg_allelic$coefficients[2,1], reg_allelic$coefficients[2,2], reg_dominant$coefficients[2,4], exp(reg_dominant$coefficients[2,1]), exp(reg_dominant$coefficients[2,1]-1.96*reg_dominant$coefficients[2,2]), exp(reg_dominant$coefficients[2,1]+1.96*reg_dominant$coefficients[2,2]), reg_dominant$coefficients[2,1], reg_dominant$coefficients[2,2])) # See header of res for elements explanation
            } # fin boucle alleles pour 1 gène
            a <- sort( as.vector(as.numeric(res[-1,6])) , decreasing = FALSE)
            p <- sort( append(p , a[1], after = length(p)) , decreasing = FALSE)  # List of Pvalue obtained for all 5 genes
          } # fin boucle pour 5 gènes
          best_Psample <- sort( append( best_Psample , p[1] , after = length(best_Psample))) ## 1000 lowest P value from each
        }
        seuil <- round(best_Psample[50] , 9) # set to 50
        seuil <- -log10(seuil)
      }



      if(!is.null(haplo) &!is.null(pcaFile) & !is.null(statut)){
        data <- utils::read.table(haplo , header= TRUE , sep = "," , dec = "." , na.strings = NA)
        statut <- stats::na.omit(utils::read.table(statut , header = TRUE , sep = "," , dec = "."))
        if(ncol(statut)!=2){stop("statut must contains 2 colunms : subjectID and affection statut.")}
        PCs <- stats::na.omit(utils::read.table(pcaFile ,  header = TRUE))
        if(ncol(pca)!=3){stop("pcaFile must contains 3 colunms : subjectID, PC1 and PC2.")}
        colnames(data)[1] <- "sample.id"
        colnames(statut)[1] <- "sample.id"
        colnames(PCs)[1] <- "sample.id"
        data <- merge(data, statut, by ="sample.id")
        all_data <- merge(data ,PCs[,1:3], by ="sample.id")
        colnames(all_data)[ncol(all_data)-2] <- "Statut"
        all_data$pheno <- ifelse(all_data$Statut=="Ctrl",0,1)
        res <- c("HLA_haplotype", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_0.95_CI_allelic", "upper_0.95_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_95%_CI_dominant", "upper_95%_CI_dominant", "Beta_dominant", "SE_dominant")

        ## Determine haplotypes to analyze; more than 2 occurences
        df <- as.data.frame(table(stats::na.omit(c( as.character(all_data$A.B.C.DR.DQ.1) , as.character(all_data$A.B.C.DR.DQ.2))))) ## number of occurences
        haplo <-subset(df , Freq >= 2) ## select haplo that appears more than 2 times
        p <- c()
        for( i in 1:1000) # set to 1000
        {
          for(hap in haplo[,1]){ # Make a loop with all the HLA haplotypes, Convert the haplotypes to 1,2 code (val1, val2), then to genotype (val3), and finally to recessive (val4) or dominant (val5) testable model
            val1 <- paste("allele1_", hap, sep="")
            all_data[,val1] <- ifelse(all_data$A.B.C.DR.DQ.1 == hap, 1,
                                      ifelse(all_data$A.B.C.DR.DQ.1 != hap, 2, NA))

            val2 <- paste("allele2_", hap, sep="")
            all_data[,val2] <- ifelse(all_data$A.B.C.DR.DQ.2 == hap, 1,
                                      ifelse(all_data$A.B.C.DR.DQ.2 != hap, 2, NA))

            val3 <- paste("genotype_", hap, sep="")
            all_data[,val3] <- paste(all_data[,val1], all_data[,val2], sep="")

            val4 <- paste("genoA_", hap, sep="")
            all_data[,val4] <- ifelse(all_data[,val3] == '11', 1,
                                      ifelse(all_data[,val3] == '12', 0,
                                             ifelse(all_data[,val3] == '21', 0,
                                                    ifelse(all_data[,val3] == '22', -1, NA))))

            val5 <- paste("genoD_", hap, sep="")
            all_data[,val5] <- ifelse(all_data[,val3] == '11', 1,
                                      ifelse(all_data[,val3] == '12', 1,
                                             ifelse(all_data[,val3] == '21', 1,
                                                    ifelse(all_data[,val3] == '22', 0, NA))))

            freq <- sum(all_data[,val1]==1,all_data[,val2]==1)/(nrow(all_data)*2) # Overall frequency
            freq_control <- sum(subset(all_data, pheno == 0)[,val1]==1,subset(all_data, pheno == 0)[,val2]==1)/(nrow(subset(all_data, pheno == 0))*2) # Controls frequency
            freq_case <- sum(subset(all_data, pheno == 1)[,val1]==1,subset(all_data, pheno == 1)[,val2]==1)/(nrow(subset(all_data, pheno == 1))*2) # Cases frequency

            reg_allelic <- summary(stats::glm(pheno ~ all_data[, val4] + PC1 + PC2 , data = all_data, family=stats::binomial(link="logit"))) # Logistic regression with allelic model
            reg_dominant <- summary(stats::glm(pheno ~ all_data[, val5] + PC1 + PC2 , data = all_data, family=stats::binomial(link="logit"))) # Logistic regression with dominant model

            res <- rbind(res, c(hap, freq, freq_control, freq_case, 4/((1/nrow(subset(all_data, pheno == 1)))+(1/nrow(subset(all_data, pheno == 0)))), reg_allelic$coefficients[2,4], exp(reg_allelic$coefficients[2,1]), exp(reg_allelic$coefficients[2,1]-1.96*reg_allelic$coefficients[2,2]), exp(reg_allelic$coefficients[2,1]+1.96*reg_allelic$coefficients[2,2]), reg_allelic$coefficients[2,1], reg_allelic$coefficients[2,2], reg_dominant$coefficients[2,4], exp(reg_dominant$coefficients[2,1]), exp(reg_dominant$coefficients[2,1]-1.96*reg_dominant$coefficients[2,2]), exp(reg_dominant$coefficients[2,1]+1.96*reg_dominant$coefficients[2,2]), reg_dominant$coefficients[2,1], reg_dominant$coefficients[2,2])) # See header of res for elements explanation

          }# fin boucle pour TOUT les haplotypes
          a <- sort( as.vector(as.numeric(res[-1,6])) , decreasing = FALSE)
          p <- sort( append(p , a[1], after = length(p)) , decreasing = FALSE)  # Autant de p value que d'haplotype ON garde une p value à chaque tour de boucle
        }
        seuil <-round(p[50] , 9) # p contains 1000 p values # set to 50
        seuil <- -log10(seuil)
      }


      if(!is.null(AA)){
        if(is.null(pcaFile)){stop(" \n pvalueCorrection function need PCA data of your cohort. ")}
        amino_acid <- utils::read.csv(AA, stringsAsFactors = F,sep=",", header = TRUE) # Import amino acid data
        amino_acid <- amino_acid[, apply(amino_acid, 2, function(x) length(unique(stats::na.omit(x)))) > 1] # Keep amino acid with information
        list_aa <- colnames(amino_acid)[-c(1:6,(ncol(amino_acid)-6):ncol(amino_acid))]

        statut <- stats::na.omit(utils::read.table(statut , header = TRUE , sep = "," , dec = "."))
        if(ncol(statut)!=2){stop("statut must contains 2 colunms : subjectID and affection statut.")}
        pca <- utils::read.table(pcaFile, header=T, stringsAsFactors = F, sep=",") # Load each HLA imputed file for each cohort
        if(ncol(pca)!=3){stop("pcaFile must contains 3 colunms : subjectID, PC1 and PC2.")}
        colnames(amino_acid)[1] <- "sample.id"
        names(statut) <- c("sample.id","Disease")
        names(pca) <- c("sample.id","PC1","PC2")

        data <- merge(amino_acid, pca[,1:3], by="sample.id")
        all_data <-  merge(data, statut, by="sample.id")
        res <- c("HLA_Amino_Acid", "Frequency", "Frequency_control", "Frequency_case", "N", "P_allelic", "OR_allelic", "lower_95%_CI_allelic", "upper_95%_CI_allelic", "Beta_allelic", "SE_allelic", "P_dominant", "OR_dominant", "lower_95%_CI_dominant", "upper_95%_CI_dominant", "Beta_dominant", "SE_dominant")
        p <- c()
        for( i in 1:1000)
        {
          print(i)
          for(aa in list_aa) # Make a loop with all amino acids
          {
            # Caculate frequency of amino acid
            freq <- mean(stats::na.omit(all_data[,aa]))/2
            freq_control <- mean(stats::na.omit(subset(all_data, Disease == 0)[,aa]))/2
            freq_case <- mean(stats::na.omit(subset(all_data, Disease == 1)[,aa]))/2

            # Convert amino acid (0,1,2) to (-1,0,1) to run the regression
            all_data[,aa][all_data[,aa] == 0] <- -1
            all_data[,aa][all_data[,aa] == 1] <- 0
            all_data[,aa][all_data[,aa] == 2] <- 1

            # Run the allelic regression
            reg_allelic <- summary(stats::glm(Disease ~ all_data[, aa] + PC1 + PC2, data = all_data, family=stats::binomial(link="logit"))) # Logistic regression with allelic model

            # Convert (-1,0,1) to (0,1) to run the regression in dominant model
            all_data[,aa][all_data[,aa] == 0] <- 1
            all_data[,aa][all_data[,aa] == -1] <- 0

            # Run the dominant regression
            reg_dominant <- summary(stats::glm(Disease ~ all_data[, aa] +  PC1 + PC2, data = all_data, family=stats::binomial(link="logit"))) # Logistic regression with allelic model
            res <- rbind(res, c(aa, freq, freq_control, freq_case, 4/((1/nrow(subset(all_data, Disease == 1)))+(1/nrow(subset(all_data, Disease == 0)))), reg_allelic$coefficients[2,4], exp(reg_allelic$coefficients[2,1]), exp(reg_allelic$coefficients[2,1]-1.96*reg_allelic$coefficients[2,2]), exp(reg_allelic$coefficients[2,1]+1.96*reg_allelic$coefficients[2,2]), reg_allelic$coefficients[2,1], reg_allelic$coefficients[2,2], reg_dominant$coefficients[2,4], exp(reg_dominant$coefficients[2,1]), exp(reg_dominant$coefficients[2,1]-1.96*reg_dominant$coefficients[2,2]), exp(reg_dominant$coefficients[2,1]+1.96*reg_dominant$coefficients[2,2]), reg_dominant$coefficients[2,1], reg_dominant$coefficients[2,2])) # See header of res for elements explanation
          }
          a <- sort( as.vector(as.numeric(res[-1,6])) , decreasing = FALSE) # Pvalues of all AA
          p <- sort( append(p , a[1], after = length(p)) , decreasing = FALSE)  # One P value for each loop

        }
        # x <- 0.05*numberLoop
        # n <- round(x,digits = -1)
        # if(n==0){n<-1}
        # seuil <- round(p[n] , 9)
        #
        seuil <- round(p[50] , 9)
        seuil <- -log10(seuil)
      }# END of AA levels




    } # END if correction == "permutation"
  } # else

  return(seuil)
}





