#!/bin/bash

# ====================================================
# 					- HLAfix -
# ====================================================
######### Script used in the R Function plink2VCF() ##########
#### Objectives:
####	* Perform QC with imputePrepSaanger pipeline provided by Sanger imputation Server
####	  tool used : plink, vcftools
####		return a VCF file ready to use in imputation server
####
####  Take 4 arguments:
#### 	1- ImputePrepSanger pipeline (done by the firstQC()function )
####  	2- Name of .bed .bim .fam without extension
#### 	3- Folder where the result will be saved
####  	4- Name of VCF file returned


if  [ $# -ne 8 ] ;
        then
      echo " Wrong number of arguments. 7 arguments are required: \n 1- imputePrepSanger.pl script 2- pathway of bed/bim/fam without extension \n 3- output folder \n 4- Name of result vcf file without extension ... \n"
      exit 1
      else
        imputePrepSanger=$1
      dataFolder=$2
      input=$3 # complete name "_postQC"
      output_folder=$4 # output_folder/imputePrep/
      filename=$5
      tabFile=$6 # HRC.r1.GRCh37.autosomes.mac5.sites.tab
      fastaFile=$7
      updateDuplicates=$8 # updateDuplicates.awk

      if [ ! -x "$output_folder" ];
      then mkdir $output_folder
      fi
      cd $output_folder/

        if [ ! -x "tmp/" ] &&  [ ! -x "clean/" ];
      then
      mkdir tmp/
        mkdir clean/
        fi


      ### ImputePrepSanger

      perl $imputePrepSanger -b $dataFolder$input".bim" -f $dataFolder$input".frq" -r $tabFile -h
      ## Run-plink.sh produced by perl but the following lines are executed instead
      mv ./Force-Allele1-$input-HRC.txt tmp/

        ### List of duplicate and removal
        plink  --bfile $dataFolder$input --write-snplist --out tmp/all_input_snps
      cat tmp/all_input_snps.snplist | sort | uniq -d > tmp/duplicated_snps.txt
      plink --bfile $dataFolder$input --exclude tmp/duplicated_snps.txt --make-bed --out clean/TEMP1

      ### Update Chromosome info and positions
      plink  --bfile clean/TEMP1 --update-chr Chromosome-$input-HRC.txt --make-bed --out clean/TEMP2
      plink  --bfile clean/TEMP2 --update-map Position-$input-HRC.txt --make-bed --out clean/TEMP3

      plink  --bfile clean/TEMP3 --flip Strand-Flip-$input-HRC.txt --make-bed --out clean/TEMP4
      plink  --bfile clean/TEMP4 --reference-allele tmp/Force-Allele1-$input-HRC.txt --make-bed --out clean/TEMP5
      plink  --bfile clean/TEMP5 --list-duplicate-vars ids-only --reference-allele tmp/Force-Allele1-$input-HRC.txt --make-bed --out clean/TEMP6

      ## Run-plink.sh
      nbrDupli=$(wc -l < clean/TEMP6.dupvar)
      if [ $nbrDupli  > 0 ];
      then
      cp clean/TEMP6.dupvar tmp/duplicatesPairs.txt
      TEMP6_dupvar=clean/TEMP6.dupvar
      if [ -f "$TEMP6_dupvar" ] && [ $(wc -l $TEMP6_dupvar) > 0 ];
      then
      plink --bfile clean/TEMP6 --extract clean/TEMP6.dupvar --reference-allele tmp/Force-Allele1-$input-HRC.txt --recode --tab --out clean/TEMP_dupliPed

      else
        plink --bfile clean/TEMP6 --reference-allele tmp/Force-Allele1-$input-HRC.txt --recode --tab --out clean/TEMP_dupliPed
      fi

      TEMP_dupliPed=clean/TEMP_dupliPed.ped
      if [[ ! -f "$TEMP_dupliPed" ]]; then
      echo "ERROR : No such file plink2VCF.sh line 69"
      exit 1
      fi
      plink --bfile clean/TEMP6 --exclude clean/TEMP6.dupvar --reference-allele tmp/Force-Allele1-$input-HRC.txt --make-bed --out clean/TEMP7


      TEMP7=clean/TEMP7.bim
      if [[ ! -f "$TEMP7" ]]; then
      echo "ERROR : No such file plink2VCF.sh line 71"
      exit 1
      fi

      awk -f $updateDuplicates clean/TEMP_dupliPed.ped > clean/TEMP_dupliPed_2.ped
      cat clean/TEMP6.dupvar | cut -d ' ' -f 2 > tmp/duplicatesRemoved.txt
      plink --file clean/TEMP_dupliPed --exclude tmp/duplicatesRemoved.txt  --reference-allele tmp/Force-Allele1-$input-HRC.txt --make-bed --out clean/TEMP8
      TEMP8=clean/TEMP8.bim
      if [[ ! -f "$TEMP8" ]]; then
      echo "ERROR : No such file plink2VCF.sh line 92"
      exit 1
      fi
      plink -bfile clean/TEMP7 --bmerge clean/TEMP8  --make-bed --out tmp/DATA_$input-updated
      DATA_input=tmp/DATA_$input-updated.bim
      if [[ ! -f "$DATA_input" ]]; then
      echo "ERROR : No such file plink2VCF.sh line 100"
      exit 1
      fi

      plink --bfile tmp/DATA_$input-updated  --a2-allele tmp/Force-Allele1-$input-HRC.txt --recode vcf bgz --keep-allele-order --out tmp/DATA_updated_vcf
      DATA_updated_vcf=tmp/DATA_updated_vcf.vcf.gz
      if [[ ! -f "$DATA_updated_vcf" ]]; then
      echo "ERROR : No such file plink2VCF.sh line 107"
      exit 1
      fi
      awk  'END {print "Number of duplicates removed: ", NR}' tmp/duplicatesRemoved.txt
      echo The duplicates removed can be found in tmp/duplicatesRemoved.txt and tmp/duplicatesPairs.txt
      else
        mv ./Force-Allele1-$input-HRC.txt tmp/
        echo Number of duplicates removed: 0 .
      fi
      fi
