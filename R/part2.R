


#######################################################
### Allele's imputation from SNPs  ###
#######################################################

#' Impute HLA allele by using HiBAG package
#' @param bedFile .bed file containing the packed binary SNP genotype data
#' @param bimFile .bim file containingthe SNP descriptions
#' @param famFile .fam file containing subject (and, possibly, family) identifiers. This is basically a tab-delimited "pedfile"
#' @param hlaLoci vector of HLA loci, by default it is: A B C DPR1
#' @param probThreshold Threshold of Postprob from Hibag to keep. By default 0.5
#' @param trainingModel Model for hibag training. Download pre-fit Model https://github.com/zhengxwen/HIBAG
#' @param pcaFile FIle with two principal components of PCA
#' @param outputFolder folder where to save the result
#' @import HIBAG
#' @importFrom utils read.table write.table
#' @importFrom stats predict
#' @importFrom snpStats read.plink
#' @importFrom stringr str_sub
#' @return All HLA alleles with imputation probabilities
#' @export
#'

hlaImputation <- function( bedFile , bimFile , famFile , hlaLoci= c("A", "B", "C", "DRB1", "DQB1"), trainingModel = NULL ,pcaFile , outputFolder, probThreshold = 0.5 )
{
  ## Check Parameters
  if (!file.exists(bedFile)){ stop("bedFile not found") }
  if (!file.exists(bimFile)){ stop("bimFile not found") }
  if (!file.exists(famFile)){ stop("famFile not found") }
  if (!file.exists(pcaFile)){ stop("pcaFile not found") }
  if (!dir.exists(outputFolder)){ stop("outputFolder directory doesn't exist.") }

  if(str_sub(bedFile,-4,-1 ) != ".bed"){ stop("bedFile doesn't have .bed extension")}
  if(str_sub(bimFile,-4,-1 ) != ".bim"){ stop("bimFile doesn't have .bim extension")}
  if(str_sub(famFile,-4,-1 ) != ".fam"){ stop("famFile doesn't have .fam extension")}

  if (is.null(trainingModel))
  { ##Build A model
     stop("Model Building not yet available. Please retry and load a model.")

  }else{
    if (!file.exists(trainingModel)) {stop("trainingModel file doesn't exist.")}
    ## Training Model preparation
    model.list <- get(load(trainingModel))

    ## HLA Prediction
    data <- HIBAG::hlaBED2Geno(bed.fn= bedFile , fam.fn = famFile, bim.fn = bimFile ) # Import Plink files

    for(hla.id in hlaLoci)
    {
      model <- HIBAG::hlaModelFromObj(model.list[[hla.id]]) # load from model, SNP of one gene (hla.id)
      pred.guess <- stats::predict(model, data, type="response", match.type="Position") # predict most probable alleles by using SNPs. SNPs of the model and dataset are matched by position.

      # Saving of best guess genotype and its postprob
      utils::write.table( pred.guess$value , paste( outputFolder , "HLA-", hla.id ,".csv",sep =""), row.names=F, quote=F, col.names=T,  sep =",",dec = "." , na = "NA")

      # # all possible combinaison of HLA alleles for each individuals
      # all_genotypes <- as.data.frame(pred.guess$postprob)
      # all_genotypes <- cbind(rownames(all_genotypes), all_genotypes)
      # colnames(all_genotypes)[1] <- "Genotypes" ## Add a colunm with combinaison of alleles  example 01:25/01:01
      # all_genotypes$allele1  <- gsub("/.*$", "",all_genotypes$Genotypes) # select allele before "/" example 01:25
      # all_genotypes$allele2  <- gsub("^.*?/", "",all_genotypes$Genotypes) # select allele after "/" example 01:01
      # df <- data.frame(all_genotypes[which(names(all_genotypes)=="allele1")] , all_genotypes[which(names(all_genotypes)=="allele2")] , all_genotypes[2: (ncol(all_genotypes)-2)] )
      # utils::write.table( df, paste( outputFolder , "HLA-", hla.id ,"_allGenotypes.csv",sep =""), row.names=F, quote=F, col.names=T,  sep =",",dec = "." , na = "NA")

    }

    plink <- snpStats::read.plink( bed = bedFile ,  bim= bimFile , fam = famFile , na.strings = "-9" )  ##return a dataframe corrresponding to the .fam
    statut <- plink$fam[,c(2,6)] # subject ID ("member") and their disease phenotype ("affected")
    names(statut) <- c("member","Disease")

    pca <- utils::read.table( pcaFile , stringsAsFactors = F, header = T, sep=',')
    names(pca) <- c("member","PC1", "PC2")
    mergedFile <- merge(statut,pca[,1:3], by="member", all = TRUE) ## add PC1, PC2 for each individuals

    for(hla.id in hlaLoci)
    {
      # Merge of the 5 files returned by hibag (best guess of genotype)
      file <- utils::read.table( paste( outputFolder , "HLA-", hla.id ,".csv",sep =""), stringsAsFactors = F, header = T, sep=',')
      file <- subset(file,  prob>=probThreshold )# keep genotype with a postprob upper than 0.5
      colnames(file) <- c("member" , hla.id ,hla.id)
      mergedFile <- merge(mergedFile , file [1:3] , by="member", all = TRUE)
    }
    names(mergedFile) <- c("subjectID" , "Disease" ,"PC1" , "PC2", "A","A","B","B","C","C","DRB1","DRB1","DQB1" ,"DQB1")
    utils::write.table( mergedFile ,paste( outputFolder, "HLAalleles_bestGuess.csv", sep="")  ,row.names=F, quote=F, col.names=T , sep = ",")
    # deletion of the 5 files HLA-A.csv , HLA-B.csv ... HLA-DQB1.csv
    system (paste("rm ", outputFolder , "HLA-*.csv ",sep =""))

}
}


#' Count HLA allele (allele dose) of interest and display their distribution across group
#'
#' Take a file containing at least colunms: subjectId , allele1 , allele 2, prob
#'
#' @param file The File containing Hla allele for each individuals and probability to have it
#' @param selection vector of allele of interest
#' @param outputFile The file that will contain the input info and
#' @return file with with all input file's data and 2xselection colunms: one for allele count(0 1 2 ) and the other with fractionned dose
#' @import ggplot2
#' @export
#'
allele_dose <- function( file , selection , outputFile )
{
  data <-  as.data.frame( utils::read.table( file , header = T, sep=',') )
  result <- as.data.frame(data$sample.id)
  colnames(result)[1] <- "sample.id"

  for(allele in  selection )
  {
    if( data$allele1 == allele )
    {
      data[, paste(allele , "_dose" , sep = "") ] <-  ifelse(data$allele2 == allele ,2,1)
      data[,paste(allele ,"_Fdose",sep = "")] <- data[,paste(allele ,"_dose",sep = "")] * data[,"prob"]
    }else{
      data[, paste( allele , "_dose" , sep = "") ] <-  ifelse(data$allele2 == allele ,1,0)
      data[,paste(allele ,"_Fdose",sep = "")] <- data[,paste(allele ,"_dose",sep = "")] * data[,"prob"]
    }
    result[,paste(allele , "_dose" , sep = "")] <-   data[, paste(allele , "_dose" , sep = "") ]
    utils::write.table(data , outputFile ,  row.names = F , sep = ",")
  }
  # Keep only Alleles Count columns And Do the total of each allele
  rownames(result) <- result[,1]
  result <- result[,-1]
  result <- rbind(result, colSums(result[,]))
  rownames(result)[nrow(result)] <- "Total.Count"

  # Replace colnames and rownames as row and columns
  result <- cbind(rownames(result),result)
  colnames(result)[1]<- "Sample.ID"
  rownames(result) <- NULL # Remove rownames

  # Subset need info for the plot
  result <- rbind(colnames(result),result)
  displ <- result[nrow(result),2:length(selection)]
  r <- as.data.frame(t(displ))
  d <- cbind(row.names(r),r)
  colnames(d) <-c("alleles" , "values")
  plot <- ggplot(d, aes(x = alleles, y = values ,fill= alleles)) + geom_bar(stat = "identity", position = position_dodge())
  print(plot)
}

