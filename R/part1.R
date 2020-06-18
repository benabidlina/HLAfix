#############################################################
###    plink2VCF : Sanger's QC and conversion into VCF    ###
#############################################################

#' Check and fix strand, check reference alleles and remove duplicate. These QC are recommended by Sanger server imputation: wrayner@well.ox.ac.uk
#'
#' @param dataFolder  folder containing bed/bim/fam and .frq files
#' @param plinkFile  bed/bim/fam file's name without extension and without pathways
#' @param outputFolder folder where to save the resul
#' @param outputFile Name wanted for your VCF
#' @importFrom snpStats read.plink col.summary
#' @importFrom R.utils countLines
#' @return VCF file ready for SNP imputation
#' @export
#'


plink2vcf <- function( dataFolder, plinkFile, outputFile, outputFolder){
  start_time <- Sys.time()
  if(!dir.exists(dataFolder)){stop("dataFolder directory not found")}
  if(!is.character(plinkFile)){stop("plinkFile must be character")}
  if(!is.character(outputFile)){stop("outputFile must be character")}

  ################ Data and scripts loading from HLAfix package  ################
  plink2VCF <- system.file("Shell", "plink2VCF.sh", package = "HLAfix") # script with bash and plink commands
  HRCcheckbim_pl <- system.file("Perl", "HRC-1000G-check-bim_v4.2.7.pl" , package="HLAfix" ) # comes with imputePrepSanger pipeline
  updateDuplicates <- system.file("awk", "updateDuplicates.awk" , package="HLAfix" ) # comes with imputePrepSanger pipeline

  ################ Download dataset
  ## fasta reference genome from 1000 genome project
  system(paste("wget -nc ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -P ", outputFolder ,sep=""))
  if(file.exists(paste(outputFolder,"human_g1k_v37.fasta.gz",sep="")))
  {
    system(paste("gunzip ", outputFolder ,"human_g1k_v37.fasta.gz",sep=""))
    fastaFile <- paste(outputFolder,"human_g1k_v37.fasta",sep="")
    system(paste("rm ", outputFolder,"human_g1k_v37.fasta.gz",sep=""))

  }else{
    stop("human_g1k_v37.fasta.gz downloading didn't work. We used the following url: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ")
  }


  ## reference genome from HRC
  system(paste("wget -nc ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1/HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz -P ", outputFolder ,sep="")) # -c for continue the download if for a reason or another it is in pause and -p for

  if(file.exists(paste(outputFolder,"HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz",sep="")))
  {
    system(paste("gunzip ", outputFolder ,"HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz",sep=""))
    tabFile <- paste(outputFolder,"HRC.r1.GRCh37.autosomes.mac5.sites.tab",sep="")
    system(paste("rm ", outputFolder,"HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz",sep=""))

  }else{
    stop("HRC.r1.GRCh37.autosomes.mac5.sites.tab downloading didn't work. We used the following url: ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1/HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz ")
  }

  ################        QC and conversion into vcf file.       ################
  system(paste("sh" , plink2VCF , HRCcheckbim_pl , dataFolder, plinkFile , outputFolder ,  outputFile , tabFile , fastaFile ,  updateDuplicates , sep = " "))
  # end_time <- Sys.time()
  # cat("\n Running time : \n")
  # end_time - start_time
}



#############################################################
###  PCA for population stratification and visualization  ###
#############################################################


#' Checking for population stratification in order to adjust and get homogenous sample
#'
#' @param bedFile .bed file containing the packed binary SNP genotype data
#' @param bimFile .bim file containingthe SNP descriptions
#' @param famFile .fam file containing subject (and, possibly, family) identifiers. This is basically a tab-delimited "pedfile"
#' @param workingFolder folder where to save the result
#' @import ggplot2
#' @importFrom grDevices dev.off  pdf
#' @importFrom utils read.table write.table
#' @return Return a plot sample's individuals ancestry with 3 references populations: European (CEU) , East-Asian (CHB) and African (YRI)
#' @export


popStrat <- function( bedFile ,  bimFile  , famFile  , workingFolder)
{

  ################ Data and scripts loading from HLAfix package  ################
  ## Reference population
  # thousgen <- system.file("extdata", "CEU_YRI_CHB.bim" , package = "HLAfix")
  thousgen <- system.file("extdata", "AFR_ASN_EUR.bim" , package = "HLAfix")
  thous <- gsub(".bim","", thousgen)
  # list_thousgen <- system.file("extdata", "ceu_yri_chb_snp.txt" , package = "HLAfix")
  list_thousgen <- system.file("extdata", "afr_asn_eur_snp.txt" , package = "HLAfix")
  system(paste("cat " , bimFile , " > " ,workingFolder ,"sample_snp_list.txt"  ,sep = "")) # Write a list of snp contained in user's file
  sample <- gsub(".bim","", bimFile) #user's file without extension


  system(paste("plink --bfile " , thous , " --extract " ,workingFolder ,"sample_snp_list.txt ", " --make-bed --out ", workingFolder , "1KG_subset" , sep = "" )) ## Extract from 1KG, snp in common with sample data
  system(paste("plink --bfile " , sample , "  --extract " , list_thousgen , " --make-bed --out ", workingFolder , "sample_subset" , sep = "" )) ## Extract from Sample, snp in common with 1kg
  system(paste("plink --bfile " , workingFolder , "sample_subset" ," --bmerge ", workingFolder , "1KG_subset"  ," --make-bed --out " , workingFolder , "merged_data" , sep = "")) ## Merge of 1KG and Sample rows

  if(file.exists(paste(workingFolder, "merged_data-merge.missnp",sep = "")))
  { # if there are SNP with 3+ alleles in the dataset, a file .missnp is created.
    system(paste("plink --bfile " , workingFolder , "sample_subset" ," --exclude ", workingFolder , "merged_data-merge.missnp"  ," --make-bed --out " , workingFolder , "sample_subset2" , sep = "")) ## Multi allelic's SNP removal
    system(paste("plink --bfile " , workingFolder , "sample_subset2" ," --bmerge ", workingFolder , "1KG_subset"  ," --make-bed --out " , workingFolder , "merged_data" , sep = "")) ## Merge of 1KG and Sample rows
    system(paste("plink --bfile " , workingFolder , "merged_data", "  --pca --out ", workingFolder, "pca_plink" , sep = ""))
    system(paste("mv ", workingFolder , "pca_plink.eigenvec" ," " , workingFolder , "pca_plink.evec" ,sep = ""))

  }else{
    system(paste("plink --bfile " , workingFolder , "merged_data", "  --pca --out ", workingFolder, "pca_plink" , sep = ""))
    system(paste("mv ", workingFolder , "pca_plink.eigenvec" ," " , workingFolder , "pca_plink.evec" ,sep = ""))
  }
  ############################  GRAPHIC REPRESENTATION OF THE ACP  ################################################
  acp <- as.data.frame(read.table(paste( workingFolder , "pca_plink.evec" , sep = ""), sep = " " , dec = ".", fill = TRUE )) # Col: FID IID followed by all different PCs
  # panelFile <- system.file("extdata", "integrated_call_samples_1kg_panel.rds", package = "HLAfix")
  # panel <- readRDS(file=panelFile)
  panelFile <- system.file("extdata", "info_pop.csv", package = "HLAfix")
  panel <- read.csv(panelFile, header = TRUE,sep = "," )

  ## check for create plot acp !!!!

  ## SAVE THE FILE THAT IS USED FOR THE PLOT
  utils::write.table(acp , paste(workingFolder , "pca_used_for_plot.csv", sep = "") , col.names = TRUE, row.names = FALSE , sep = "," )
  cat(paste(workingFolder , "pca_used_for_plot.csv is saved. \n ", sep = ""))

  message("Warning: \n Ancestry identification of individuals is done with plink but it is less precise than EigenSoft and doesn't remove outliers automatically. \n If you want you can use Eigensoft and be sure of your sample's homogenesity before continuing \n")
  message("Please remove outliers and submit bed bm fam files with homogeneity ")

  # ANCESTRY PLOT
  # acp_plot<- ggplot2::ggplot(selection , aes(x=PC1, y=PC2, color=group)) + geom_point() + ggtitle("Ancestry plot ") + xlab("PC1") + ylab("PC2")
  data_1kg <- selection[(selection[,"group"]=="EUR" |selection[,"group"]=="East-Asian" | selection[,"group"]=="AFR"| selection[,"group"]=="Indian"),]
  data_sample <- selection[selection[,"group"]=="sample",]
  # Counts
  nSample <- nrow(data_sample)
  nAFR <- nrow(selection[selection[,"group"]=="AFR",])
  nASN <- nrow(selection[selection[,"group"]=="East-Asian",])
  nEUR <-  nrow(selection[selection[,"group"]=="EUR",])
  nIndian <- nrow(selection[selection[,"group"]=="Indian",])

  acp_plot <- ggplot2::ggplot(data_sample , aes(x=PC1, y=PC2, color=group)) + geom_point(data=data_1kg , aes(x=PC1,y=PC2,color=group) , alpha = 0.2,size=2) + geom_point(size=2) + ggtitle("Ancestry plot ") + xlab("PC1") + ylab("PC2") +  theme(plot.title = element_text(hjust = 0.5)) +  scale_color_manual(labels = c(paste("AFR (N=",nAFR, ")",sep = ""), paste("East-Asian (N=",nASN, ")",sep = "") , paste("EUR (N=",nEUR, ")",sep = ""), paste("Indian (N=",nIndian, ")",sep = ""), paste("Sample (N=",nSample, ")",sep = "")) ,values = c("#f3ab55","#65da11","#1563ea","#f244da","#c682ec" ))

  pcaFile <- acp[,c(2,3,4)] # col IID, PC1 and PC2
  pcaFile[,1] <- gsub("NA\\d{5}","1KG",pcaFile[,1])
  pcaFile[,1] <- gsub("HG\\d{5}","1KG",pcaFile[,1])
  pcaFile <- pcaFile[pcaFile[,1]!="1KG",] # exclude info concerning 1KG individuals
  colnames(pcaFile)<- c("sample.id","PC1","PC2")
  utils::write.table(pcaFile , paste(workingFolder , "pcaFile.csv", sep = "") , col.names = TRUE, row.names = FALSE , sep = "," )
  grDevices::pdf(file = paste(workingFolder , "Ancestry plot.pdf", sep = ""))
  print(acp_plot)
  grDevices::dev.off()
  rm(panel) # remove from env the panel data
  # remove intermediary files created
  system(paste("rm ",workingFolder,"1KG_subset* ",workingFolder,"merged_data* ",workingFolder,"sample_subset* ",workingFolder,"sample_snp_list.txt " , sep = "" ))
}
