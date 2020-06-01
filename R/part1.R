if(getRversion() >= "2.15.1") utils::globalVariables(".")
globalVariables(c("V3","V4" , "group","too_close","values","green","alleles","PC1","PC2"))

#######################################################
###  QC for imputation & population stratification  ###
#######################################################

#' Do a standardised quality control, it cleans your data by removing not informative's SNP and analyze stratification of your population
#' @param bedFile .bed file containing the packed binary SNP genotype data
#' @param bimFile .bim file containingthe SNP descriptions
#' @param famFile .fam file containing subject (and, possibly, family) identifiers. This is basically a tab-delimited "pedfile"
#' @param vcf variant call format file
#' @param outputFolder folder where to save the result
#' @param outputFile Name wanted for your VCF
#' @param shFile absolute path of the bash script for quality control
#' @importFrom snpStats read.plink col.summary
#' @importFrom R.utils countLines
#' @return Write a summary of the data after the quality control and cleaned file
#' @export
#'

firstQC <- function( bedFile=NULL , bimFile=NULL , famFile=NULL ,vcf=NULL, outputFolder , outputFile , shFile = NULL)
{
  # Check Parameters
  if(!file.exists(bedFile)){stop("bedFile is not found ")}
  if(!file.exists(bimFile)){stop("bimFile is not found ")}
  if(!file.exists(famFile)){stop("famFile is not found ")}
  if(!is.character(outputFile)){stop("outputFile must be character")}

  if(is.null(bedFile) & is.null(bimFile) & is.null(famFile) & is.null(vcf)){ stop("Data not provided. Please provide either bed,bim,fam or vcf to proceed.")}
  if(is.null(bedFile) & !is.null(bimFile) & !is.null(famFile) & is.null(vcf)){ stop("One plink file is missing : bedFile not given.")}
  if(!is.null(bedFile) & is.null(bimFile) & !is.null(famFile) & is.null(vcf)){ stop("One plink file is missing : bimFile not given.")}
  if(!is.null(bedFile) & !is.null(bimFile) & is.null(famFile) & is.null(vcf)){ stop("One plink file is missing : famFile not given.")}
  if(is.null(bedFile) & is.null(bimFile) & !is.null(famFile) & is.null(vcf)){ stop("Two plink file is missing : bedFile and bimFile not given.")}
  if(is.null(bedFile) & !is.null(bimFile) & is.null(famFile) & is.null(vcf)){ stop("Two plink file is missing : bedFile and famFile not given.")}
  if(!is.null(bedFile) & is.null(bimFile) & is.null(famFile) & is.null(vcf)){ stop("Two plink file is missing : bimFile and famFile not given.")}

  if(!is.null(bedFile) & !is.null(bimFile) & !is.null(famFile) & is.null(vcf)){
    if(!file.exists(bedFile)){ stop("File passed in bedFile parameter doesn't exist ")}
    if(!file.exists(bimFile)){ stop("File passed in bimFile parameter doesn't exist ")}
    if(!file.exists(famFile)){ stop("File passed in famFile parameter doesn't exist ")}
  }
  if(!is.character(outputFile)){ stop("outputFile must be string ")}
  if(!is.null(vcf)){
    system("plink --vcf ",vcf," --make-bed --out ", "inputData")
    system(paste("plink --vcf ", vcf , " --make-bed --out ", outputFolder, "firstInput" , sep = ""))
    bedFile <- paste(outputFolder,"firstInput.bed",sep = "")
    bimFile <- paste(outputFolder,"firstInput.bim",sep = "")
    famFile <- paste(outputFolder,"firstInput.fam",sep = "")
  }


  ################ Quality control execution via bash script ################
  mind <- readline(prompt = " \n SNP missingness by individual (mind)? ")
  if( mind < 0 || mind > 1 ) {stop("mind must be a value between 0-1")}
  geno <- readline(prompt = "\n SNP missingness in genotype (geno)? ")
  if( geno < 0 || geno > 1 ) {stop("geno must be a value between 0-1")}
  maf <- readline(prompt = "\n SNP minor allelic frequency threshold (maf)? ")
  if( maf < 0 || maf > 1 ) {stop("maf must be a value between 0-1")}
  hwe <- readline(prompt = "\n SNP Hardy weinberg equilibrium threshold (hwe)? ")
  if( hwe < 0 || hwe > 1 ) {stop("hwe must be a value between 0-1")}

  system2( command =  shFile ,args = c( bimFile , outputFolder ,outputFile, mind, geno , maf , hwe ))

  ################ Genotype Summary Before QC  ################
  plinkFile <- gsub(".bim","", paste(outputFolder,"/tmp/input_goodID_NOduplicate.bim")) #user's file without extension
  if(file.exists(paste(outputFolder,"/tmp/input_goodID_NOduplicate.bim",sep = "")))
  {
    genotype <- genio::read_plink(paste( outputFolder , "/tmp/input_goodID_NOduplicate" , sep = "") )
  }else{
    genotype <- genio::read.plink(paste(outputFolder,"/tmp/input_no_duplicate" ,sep = "") )
  }

  statutDisease <- genotype$fam[,c(2,6)]
  ctrlNumber <-  nrow(statutDisease[statutDisease[,2]==1,2])
  if(is.null(ctrlNumber)){ctrlNumber <- 0}
  caseNumber <-  nrow(statutDisease[statutDisease[,2]==2,2])
  if(is.null(caseNumber)){caseNumber <- 0}
  indNumber <- R.utils::countLines(paste(outputFolder,"/tmp/input_no_duplicate.fam" ,sep = ""))
  snpNumber <- R.utils::countLines(paste(outputFolder,"/tmp/input_no_duplicate.bim" ,sep = ""))

  # ### Small table
  cat( "\n",
       " .-------------------------------------------.\n",
       "|            Summary of input data           |\n",
       "+-------------------------+------------------+\n"  ,
       "|  Number of individuals  |  Numbers of SNP  |\n",
       "+-------------------------+------------------+\n"  ,
       "|         ", indNumber , "                               |\n",
       "+------------+------------+------------------+\n" ,
       "|  Control   |     Case   |                  |\n",
       "+------------+------------+------------------+\n" ,
       "|    ", ctrlNumber , "         ", caseNumber ,"             ", snpNumber , "    |\n",
       "+-------------------------+------------------+\n \n"   )



  ################ Qualty Control and summary of big steps  ################
  if(isTRUE(file.exists(paste( outputFolder, "/clean/clean1.bim" , sep = ""))) && isTRUE(file.exists(paste( outputFolder, "/clean/clean2.bim" , sep = ""))) && isTRUE(file.exists(paste( outputFile, "_postQC.bim" , sep = "")))){
    clean1_snp <- R.utils::countLines(paste( outputFolder, "/clean/clean1.bim" , sep = ""))
    clean1_snp <- gsub("[A-Z]/.,:", "", clean1_snp)
    clean1_ind <- R.utils::countLines(paste( outputFolder, "/clean/clean1.fam" , sep = ""))
    clean1_ind <- gsub("[A-Z]/.,:", "", clean1_ind)
    clean2_snp <- R.utils::countLines(paste( outputFolder, "/clean/clean2.bim" , sep = ""))
    clean2_snp <- gsub("[A-Z]/.,:", "", clean2_snp)
    clean2_ind <- R.utils::countLines(paste( outputFolder, "/clean/clean2.fam" , sep = ""))
    clean2_ind <- gsub("[A-Z]/.,:", "", clean2_ind)
    postQC_snp <- R.utils::countLines(paste(outputFile , "_postQC.bim" , sep = ""))
    postQC_snp <- gsub("[A-Z]/.,:", "", postQC_snp)
    postQC_ind <- R.utils::countLines(paste(outputFile , "_postQC.fam" , sep = ""))
    postQC_ind <- gsub("[A-Z]/.,:", "", postQC_ind)
    if(file.exists(paste( outputFolder, "tmp/snp_chrXYMT_to_exlude.txt" , sep = "")))
    {
      snp_chr <- R.utils::countLines(paste( outputFolder, "tmp/snp_chrXYMT_to_exlude.txt" , sep = ""))
      snp_chr <- gsub("[A-Z]/.,:", "", snp_chr)

      if( snp_chr > 0 ){cat( "\n " , snp_chr  , "SNPs on chromosome X, Y , mitochondrial and/or bi-allelic SNPs removed \n")}

    }
    if(file.exists(paste( outputFolder, "tmp/too_related_id.txt", sep = "")))
    {
      too_close <- R.utils::countLines(paste(outputFolder, "tmp/too_related_id.txt" , sep = ""))
      who <- read.table(paste(outputFolder, "tmp/too_related_id.txt" , sep = ""), stringsAsFactors=FALSE)
      l <- who$V2
      cat(" \n ", too_close , " individual(s) removed after --genome : \n")
      for(i in 1:length(who$V2)){ print(l[i])};
      cat("\n" )
    }

    cat("--mind ", mind , " --geno " , geno ," --hwe ", hwe  ," --maf ", maf ,
        "\n ", postQC_snp ," SNPs and " , postQC_ind , " peoples remained. \n")

  }

  ################ Genotype Summary after QC  ################
  if(file.exists(paste(outputFile , "_postQC.bed" , sep = "")))
  {
    geno <- genio::read_plink(paste( outputFile , "_postQC" , sep = "") )
    snp <- R.utils::countLines(paste( outputFolder, "/clean/clean1.bim" , sep = ""))
    snp <- gsub("[A-Z]/.,:", "", clean1_snp)
    ind <-  R.utils::countLines(paste(outputFolder, "/clean/clean1.fam" , sep = ""))
    ind <- gsub("[A-Z]/.,:", "", clean1_ind)

    cat( "\n",
         "@--------------------------------------------@\n",
         "|            Summary of post QC data         |\n",
         "+-------------------------+------------------+\n"  ,
         "|  Number of individuals  |  Numbers of SNP  |\n",
         "+-------------------------+------------------+\n"  ,
         "|         ", indNumber , "                    " , snp ,"   |\n",
         "@-------------------------+------------------@\n \n"   )

  }
}



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

  # Create a colunm with pop name for 1KG individuals
  for(x in panel[!is.na(panel[,1]),1] )
  { x <- as.character(x)
  acp[acp[,2]==x,23] <- panel[panel[,1]==x,3]
  }

  # Recode ASN, YRI, EUR into CHB YRI and CEU
  acp[,1] <- ifelse(acp$V23=="East_ASN", "East-Asian",
                    ifelse(acp$V23=="South_ASN", "Indian",
                           ifelse(acp$V23=="AFR", "AFR",
                                  ifelse(acp$V23=="EUR", "EUR","sample")))) # remplace col FIID by the group

  # Recode what left into "sample" group
  suppressWarnings({
    acp[is.na(acp[,1]),1] <- "sample" # it's working !
  }) # "sample" is well written for individuals coming from the sample


  # Warning message
  names(acp)[1] <- "group"
  names(acp)[2] <- "IID"
  names(acp)[3] <- "PC1"
  names(acp)[4] <- "PC2"
  selection <- acp[,c("group","PC1","PC2")]
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





















