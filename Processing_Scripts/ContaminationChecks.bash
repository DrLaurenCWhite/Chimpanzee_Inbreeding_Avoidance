#!/bin/bash


#Before merging libraries across sequencing runs: check for contamination using Kraken and Heterozygosity estimates
#Before merging libraries from the same individuals: Check for sample or ID mix-ups by checking relatedness. Libraries from the same individual should have R=1


##Build Kraken database
#Requires a reference genomes in fasta format
mkdir ./KrakenContamDB
kraken2-build --download-taxonomy --db ./KrakenContamDB
kraken2-build --download-library human --db ./KrakenContamDB
for dir in ./FastaGenome/Primates/*/
do
	if [[ ${dir} != *"Human"* ]] #dont need the human reference as we downloaded it from ncbi above
	then
		echo "Add to library"
		kraken2-build --add-to-library ${dir}* --db ./KrakenContamDB
	fi
done
kraken2-build --build --threads 10 --db ./KrakenContamDB
kraken2-build --clean --db ./KrakenContamDB


##Run Kraken for each library (Before merging)
#Requires mapped reads for each library in fasta format
for i in ./FastaFiles/*_pT6_mapped.fa
do
	samp=`echo ${i} | cut -d"_" -f1 | cut -d"/" -f3`
	echo ${samp}
	kraken2 --threads 5 --db ./KrakenContamDB/ --report ./Reports/${samp}_report.txt  --report-zero-counts --use-names ./FastaFiles/${samp}_pT6.fa > ./Output/${samp}_output.txt
done


##Get heterozygosity estimates per library (Before merging)
depth=3
for i in *_pT6_mapped.bam
do
	echo ${i}
	samp=`echo ${i} | cut -d"_" -f1`
	echo ${samp}
	angsd -P 1 -i ./${i} -rf /mnt/scratch/lauren/Contam/Chromosomes.txt -anc /mnt/scratch/lauren/FastaGenome/panTro5.fa -ref /mnt/scratch/lauren/FastaGenome/panTro5.fa -out /mnt/scratch/lauren/Contam/SFS_Files/${samp}_md${depth} -dosaf 1 -fold 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepth ${depth} -setMaxDepth 200 -doCounts 1 -GL 1 
	realSFS ${samp}_md${depth}.saf.idx -P 6 > ${samp}_md${depth}.sfs 
	#in R: extract het
	#a<-scan("${samp}_md${depth}.sfs")
	#a[2]/sum(a)
done



##ngsRelate
#Just run on chr1 for efficiency
#Requires a list of libraries to merge (list excludes libraries identified as potentially contaminated from above steps)
minind=50
maf=0.05
MAF="05"
CHROM="chr1"
angsd -P 1 -b ./AllLibs_DeContam.bamlist -r ${CHROM} -ref /panTro6.fa -out ./IDandSexCheck/AllLibs_DeContam_minInd${minind}_maf${MAF}_${CHROM} -gl 2 -doGlf 3 -doMajorMinor 1 -doMaf 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -SNP_pval 1e-6 -minInd ${minind} -minmaf ${maf} 
zcat ./IDandSexCheck/AllLibs_DeContam_minInd${minind}_maf${MAF}_${CHROM}.mafs.gz | cut -f6 | sed 1d > ./IDandSexCheck/AllLibs_DeContam_minInd${minind}_maf${MAF}_${CHROM}.freq
ngsRelate -g ./IDandSexCheck/AllLibs_DeContam_minInd${minind}_maf${MAF}_${CHROM}.glf.gz -n 783  -f ./IDandSexCheck/AllLibs_DeContam_minInd${minind}_maf${MAF}_${CHROM}.freq  -O ./IDandSexCheck/AllLibs_DeContam_minInd${minind}_maf${MAF}_${CHROM}

