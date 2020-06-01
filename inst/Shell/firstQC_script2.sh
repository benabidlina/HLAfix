#!/bin/bash
# ====================================================
# 					- HLAfix -
# ====================================================
######### Script used in the R Function firstQC() ##########
#### Description:
####	 Perform quality control using plink analysis tool
####		return 3 files .bed .bim .fam for population stratification visualization popStrat()
#### 		return 3 files .bed .bim .fam cleaned and ready for plink2vcf() function
####
#### Take 4 arguments:
#### 	1- Pathway of bim file (bed/bim/fam must have the same name and in the same folder)
#### 	2- Folder where the result will be saved
####  	3- Name of result file without extension
####  	4- global missingness per individuals (--mind)
####  	5- missingness in genotype (--geno)
####  	6- Minor allelic frequency (--maf)
####  	7- Hardy weinberg equilibrium threshold (--hwe)

plink(){

/usr/bin/plink1.9 "$@"
}
export plink

x=$1
y=$(echo $x | tail -c 5)

if  [ $# -ne 7 ] ;
	then echo " Wrong number of arguments. 3 arguments are required: \n 1- pathway of bim file (bed/bim/fam must have the same) \n 2- output folder \n 3- Name of result file without extension \n"
	exit 1

else
	#echo "${x: -4}" != ".bim"   #check what are the 4 last character (the extension)

	if [ $y != ".bim"  ];
		then echo " First agument's extension forgotten or incorrect"
		exit 1
	else
		input=$(echo "$x" | cut -d '.' -f 1)
		output_folder=$2
		filename=$3
		mind=$4
		geno=$5
		maf=$6
		hwe=$7


		if [ -f "$input" ];
			then echo "FILE not found"
			exit 1
		else
				if [ ! -x "$output_folder" ];
					then mkdir $output_folder

				fi

				cd $output_folder/

				if [ ! -x "tmp/" ] &&  [ ! -x "clean/" ];
				then
					mkdir tmp/
					mkdir clean/
				fi

				chrSexNbr=$(awk '{ if($1=="23" || $1=="24" || $1=="25 || $1=="26) print $1 " " $2 " " $3 ;}' "$input.bim" | wc -l)
        noRSID=$(awk '{ if($2==".") print $1 " " $2 " " $3 ;}' "$input.bim" | wc -l)

				if [ noRSID > 0 ] ;
				then
				 ## Set missing variants
					plink --bfile $input --set-missing-var-ids '@:#' --make-bed --out tmp/input_goodID
				fi
				 ## Gather SNP with A/T or C/G alleles and SNP with indels
				awk '($5=="A" && $6=="T")||($5=="T" && $6=="A")||($5=="C" && $6=="G")||($5=="G" && $6=="C") {print $2}' "tmp/input_goodID.bim" > tmp/list_AT_GC.txt
			 	cat "tmp/input_goodID.bim" | awk '($5=="-" || $6=="-" ) {print $2}'  > tmp/indel_to_exclude.txt

				if [ chrSexNbr > 0 ] ;
				then
					## Gather SNP on Chr XYMito
				grep -w 23 "$input.bim" | awk '{print $2}' > tmp/chsX.txt
				grep -w 24 "$input.bim" | awk '{print $2}' > tmp/chs24.txt
				grep -w 25 "$input.bim" | awk '{print $2}' > tmp/chs25.txt
				grep -w 26 "$input.bim" | awk '{print $2}' > tmp/chs26.txt
				cat tmp/chsX.txt tmp/chs24.txt tmp/chs25.txt tmp/chs26.txt tmp/list_AT_GC.txt > tmp/snp_chrXYMT_to_exlude.txt

					## Check sex
					if [ $(awk '{ if($1=="23") print $2;}' "$input.bim" | wc -l) > 0 ] ;
					then
					plink --bfile $input --check-sex --out tmp/sex_check
					fi
				fi

			 	## Remove duplicates for summary and next steps for QC
			  uniq -u "tmp/input_goodID.bim" | uniq --unique --check-chars=18  > "tmp/duplicate_check.bim"
			  diff "tmp/input_goodID.bim"  "tmp/duplicate_check.bim" > "tmp/duplicate_to_exclude.txt"
			  plink --bfile tmp/input_goodID --exclude tmp/duplicate_to_exclude.txt --make-bed --out tmp/input_goodID_NOduplicate

        uniq -u $input".bim" | uniq --unique --check-chars=18  >  "tmp/duplicate_check_input.bim"
			  diff $input".bim" "tmp/duplicate_check_input.bim" > "tmp/duplicate_to_exclude_bim.txt"
			  plink --bfile $input --exclude tmp/duplicate_to_exclude_bim.txt --make-bed --out tmp/input_no_duplicate

			  ## INDELS deletions
				plink --bfile tmp/input_goodID_NOduplicate --exclude tmp/indel_to_exclude.txt --make-bed --out $input"_noINDEL"
				plink --bfile tmp/input_goodID_NOduplicate --exclude tmp/snp_chrXYMT_to_exlude.txt --make-bed --out clean/clean1

CLEAN1=clean/clean1.bim
if [[ ! -f "$CLEAN1" ]]; then
	echo "ERROR : Step 1 of QC failed to exclude SNPs on chromosomes : X, Y and MT , and to remove SNPs duplicates ."
	exit 1
fi

				 ## Exclude SNPs with bad QC
				plink --bfile clean/clean1 --mind $mind --geno $geno --maf $maf --hwe $hwe --make-bed --out clean/clean2

CLEAN2=clean/clean2.bim
if [[ ! -f "$CLEAN2" ]]; then
	echo "ERROR : Step 2 of QC failed to exclude SNPs having bad quality --mind 0.05 --geno 0.02 --maf 0.01 --hwe 0.001"
	exit 1
fi

				plink --bfile clean/clean2 --genome --out tmp/relationship
				cat tmp/relationship.genome | awk '{if ($10 > 0.125) {print $1 " " $2;}}'| tail -n 2 > tmp/too_related_id.txt

				 ## 3- remove too related ind listed in too_related_id.txt

				plink --bfile clean/clean2 --remove  tmp/too_related_id.txt --make-bed --out  $filename"_postQC"

				# map & ped files

			 	plink --bfile $filename"_postQC" --recode --out $filename"_postQC"

			 	# SNP frequencies

			 	plink --bfile $filename"_postQC" --freq --out $filename"_postQC"

FILE_result="$filename"
Ext_FILE_result="_postQC.frq"
FILE_result="$FILE_result$Ext_FILE_result"
if [[ ! -f "$FILE_result" ]]; then
	echo "ERROR : Step 3.1 of QC failed to calculate MAF frequency & to recode files into ped and map & to remove too related individuals "
	exit 1
fi

				 ## 4- QC for stratification

					# >Exclude SNPs with r2>0.2
					# >Give an output filename with SNP to exclude
				plink --bfile $filename"_postQC" --indep-pairwise 50 5 0.2
				plink --bfile $filename"_postQC" --exclude plink.prune.out --make-bed --out clean/clean3
CLEAN3=clean/clean3.bim
if [[ ! -f "$CLEAN3" ]]; then
	echo "ERROR : Step 3.2 of QC failed to exclude SNPs with r2>0.2"
	exit 1
fi
## 5- Select SNPs in a range (MHC region)
				plink --bfile clean/clean3 --chr 6 --from-mb 29 --to-mb 34 --write-snplist --out tmp/plink
					# *Output is in plink.snplist
				plink --bfile clean/clean3 --exclude tmp/plink.snplist --make-bed --out clean/clean4
CLEAN4=clean/clean4.bim
if [[ ! -f "$CLEAN4" ]]; then
	echo "ERROR : Step 4 of QC failed to extract SNPs in MHC region"
	exit 1
fi

				# 6- Use maf 0.1
				plink --bfile clean/clean4 --maf 0.1 --make-bed --out $filename"_postStratQC"
FILE_result_CLEAN5="$filename"
Ext_FILE_result_CLEAN5="_postStratQC.bim"
CLEAN5="$FILE_result_CLEAN5$Ext_FILE_result_CLEAN5"
if [[ ! -f "$CLEAN5" ]]; then
	echo "ERROR : Step 5 of QC failed to finalize pop stratification QC ( --maf 0.1 )"
	exit 1
fi
				mv plink* tmp/
				# rm -R tmp/ clean/


		fi
	fi
fi

