###connect to computer MW
ssh dpd64@128.84.15.23
cd /Volumes/ManihotDB/June2016_VCF
###send the vcf file to cbsu
scp -r /Volumes/ManihotDB/June2016_VCF/RawVCF rjl278@128.84.3.220:/home/rjl278
###send names clones Lydia to cbsu
scp -r ~/Google\ Drive/GWAS_lydia/clonesGBS.txt rjl278@128.84.3.220:/home/rjl278
###cbsu
ssh -X rjl278@cbsumm07.tc.cornell.edu
cd /workdir
cp -r /home/rjl278/RawVCF RawVCF
cp -r /home/rjl278/clonesGBS.txt clonesGBS.txt
# Write June2016 Raw VCF sample list to file
for i in {1..3}; do vcftools --gzvcf RawVCF/cassavaGBSbuild_June2016_withRef_chr${i}.vcf.gz --keep clonesGBS.txt --recode --recode-INFO-all --stdout | bgzip -c > Lydia_build_chr${i}.vcf.gz ; done
for i in {4..6}; do vcftools --gzvcf RawVCF/cassavaGBSbuild_June2016_withRef_chr${i}.vcf.gz --keep clonesGBS.txt --recode --recode-INFO-all --stdout | bgzip -c > Lydia_build_chr${i}.vcf.gz ; done
for i in {7..9}; do vcftools --gzvcf RawVCF/cassavaGBSbuild_June2016_withRef_chr${i}.vcf.gz --keep clonesGBS.txt --recode --recode-INFO-all --stdout | bgzip -c > Lydia_build_chr${i}.vcf.gz ; done

######################################################################################
#transfer file to rjl computer
######################################################################################

scp -r ~/Google\ Drive/GWAS_lydia/clonesGBS.txt dunia@128.253.192.28:/home/roberto/Desktop/Lydia
cp -r /home/DB2/GBS_June_2016 GBS_June_2016
cp -r home/DB2/dunia/clonesGBS.txt clonesGBS.txt
cp -r /home/roberto/Software/beagle.05Jul16.587.jar beagle.05Jul16.587.jar

#test do only once check number of individuals take any chromosome vcf
vcftools --gzvcf  GBS_June_2016/chr4.vcf.gz --keep clonesGBS.txt --recode --stdout | gzip -c > Lydia_chr4.vcf.gz
gunzip Lydia_chr4.vcf.gz
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' Lydia_chr4.vcf

#subset vcf file
for i in {1..18}; do vcftools --gzvcf GBS_June_2016/chr${i}.vcf.gz --keep clonesGBS.txt --recode --stdout | gzip -c > Lydia_chr${i}.vcf.gz; done
cd /home/roberto/Desktop/Lydia
#concatenate subsets of vcf files 
vcf-concat Lydia_chr1.vcf.gz Lydia_chr2.vcf.gz Lydia_chr3.vcf.gz Lydia_chr4.vcf.gz Lydia_chr5.vcf.gz Lydia_chr6.vcf.gz Lydia_chr7.vcf.gz Lydia_chr8.vcf.gz Lydia_chr9.vcf.gz Lydia_chr10.vcf.gz Lydia_chr11.vcf.gz Lydia_chr12.vcf.gz Lydia_chr13.vcf.gz Lydia_chr14.vcf.gz Lydia_chr15.vcf.gz Lydia_chr16.vcf.gz Lydia_chr17.vcf.gz Lydia_chr18.vcf.gz | gzip -c > Lydia.vcf.gz
# Filter 1
vcftools --gzvcf Lydia.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | awk '$4 != "-" {print}' | awk '$5 != "-" {print}' | gzip -c > Lydia.vcf.gz_filter1.vcf.gz; 
# Filter 2
vcftools --gzvcf Lydia.vcf.gz_filter1.vcf.gz --minDP 1 --maxDP 50 --recode --recode-INFO-all --stdout | gzip -c > Lydia.vcf.gz_filter2.vcf.gz
# Filter 3
vcftools --gzvcf Lydia.vcf.gz_filter2.vcf.gz --remove-filtered-all --max-missing 0.4 --max-meanDP 120 --thin 5 --recode --stdout | gzip -c > Lydia.vcf.gz_filter3.vcf.gz
#running the imputation

#running in CBSU
wget http://faculty.washington.edu/browning/beagle/beagle.23Jul16.fb0.jar
export JAVA_HOME=/usr/local/jdk1.8.0_45; export PATH=$JAVA_HOME/bin:$PATH; 
java -Xms2g -Xmx50g -jar beagle.23Jul16.fb0.jar gl=Lydia.vcf.gz_filter3.vcf.gz out=Imputed.Lydia.vcf.gz_filter3.vcf.gz nthreads=40 window=3200 overlap=320 niterations=10 
#running in RJL
java -Xms2g -Xmx50g -jar beagle.05Jul16.587.jar gl=Lydia.vcf.gz_filter3.vcf.gz out=Imputed.Lydia.vcf.gz_filter3.vcf.gz nthreads=12 window=3200 overlap=320 niterations=10 
gunzip Lydia.vcf.gz_filter3.vcf.gz
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' Lydia.vcf.gz_filter3.vcf
cat Imputed.Lydia.vcf.gz_filter3.vcf | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c #214759 TYPE=snp
gzip Imputed.Lydia.vcf.gz_filter3.vcf

scp -r rjl278@cbsulm04.tc.cornell.edu:/workdir/Imputed.Lydia.vcf.gz_filter3.vcf.gz.vcf.gz ~/Google\ Drive/GWAS_lydia/Imputed.Lydia.vcf.gz_filter3.vcf.gz.vcf.gz

####FILTER BY AR needs to be bgzipped first
###import Imputed.Lydia.vcf.gz_filter3.vcf.gz
scp -r dunia@128.253.192.28:/home/roberto/Desktop/Lydia/Imputed.Lydia.vcf.gz_filter3.vcf.gz   ~/Google\ Drive/GWAS_lydia
gunzip Imputed.Lydia.vcf.gz_filter3.vcf.gz
/Users/dpd64/htslib/bgzip Imputed.Lydia.vcf.gz_filter3.vcf
/Users/dpd64/htslib/tabix -p vcf Imputed.Lydia.vcf.gz_filter3.vcf.gz && /Users/dpd64/bcftools/bcftools view --include 'INFO/AR2>=0.3' Imputed.Lydia.vcf.gz_filter3.vcf.gz | gzip -c > Imputed.Lydia._AR2filtered.vcf.gz 

# How many SNPs (at MAF > 1%)
vcftools --gzvcf Imputed.Lydia._AR2filtered.vcf.gz  --recode --maf 0.01 --stdout | gzip -c > Imputed.Lydia._AR2filtered_maf.vcf.gz

#write a plink file
gunzip Imputed.Lydia._AR2filtered_maf.vcf.gz
./plink --vcf Imputed.Lydia._AR2filtered_maf.vcf --out Imputed.Lydia._AR2filtered_maf
./plink --bfile Imputed.Lydia._AR2filtered_maf --recode --tab --out Imputed.Lydia._AR2filtered_maf

# Write a dosage file for further analysis  
vcftools --vcf Imputed.Lydia._AR2filtered_maf.vcf --extract-FORMAT-info DS --out Dosage.Imputed.AR.maf.Lydia.vcf

#Import
scp -r dunia@128.253.192.28:/home/roberto/Desktop/Lydia/Imputed.Lydia.vcf.gz_filter3.vcf.gz   ~/Google\ Drive/GWAS_lydia

# Load in R
dose<-read.table("Dosage.Imputed.AR.maf.Lydia.vcf", header=T, stringsAsFactors = F)
dim(dose)
ids<-paste(paste0("S",dose$CHROM),dose$POS,sep="_")
dose1<-t(dose[,-c(1,2)])
colnames(dose1)<-ids
rm(dose); gc()
snps<-dose1;
save(snps,file="LydiaSamples_AllChrom.Rdata")
           
