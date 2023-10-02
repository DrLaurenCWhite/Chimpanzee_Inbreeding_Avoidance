#!/bin/bash

##Processing steps used to generated relatedness estimates

##Demultiplex each file (already seperated by index (i.e. read group)###
#Requires barcode file in the form of <barcode> <sample_R1.fa> <sample_R2.fa>
Run="Sequencing Run Name"
RG=`list of read group names`
for i in ${RG}; do #Demux on P7 barcode
	echo ${i}
	cd ./${i}
	mkdir -p ./2_demux/
	rm ./2_demux/*
	echo "Demultiplexing on BC2"
	sabre pe -m 1 -f <(unpigz -c ./1_raw/*R2*.fastq.gz) -r <(unpigz -c ./1_raw/*R1*.fastq.gz) -b ./${Run}/barcodes/${i}_sabre_BC2.txt -u ${i}_unknown_demux1_R2.fastq -w ${i}_unknown_demux1_R1.fastq > ./2_demux/${i}_Report_demux1.txt
	mv *fastq ./2_demux/
	cd ../
done
for i in ${RG}; do #Demux on P5
	echo ${i}
	cd ./${i}
	echo "Demultiplexing on BC1"
	if [[ -f ./2_demux/${i}_Report_demux2.txt ]]; then
		rm ./2_demux/${i}_Report_demux2.txt
	fi
	for READS1 in ./2_demux/*_R1_demux1.fastq; do
		SAMP=`echo $(basename "$READS1") | cut -d "_" -f 1`
		echo $SAMP
		sabre pe -m 1 -f ${READS1} -r ${READS1/R1/R2} \
		-b ./${Run}/barcodes/${i}/${SAMP}_sabre_BC1.txt -u ${SAMP}_unknown_demux2_R1.fastq -w ${SAMP}_unknown_demux2_R2.fastq >> ./2_demux/${i}_Report_demux2.txt
		mv *fastq ./2_demux/
	done
	pigz ./2_demux/*fastq
	cd ../
done


##Adapter trimming
#Requires fasta file of adapter sequences and a list of all samples per read group
Run="Sequencing Run Name"
RG=`list of read group names`
bbduk=path to bbduk.sh
for i in ${RG}; do
	while read SAMPLE; do
		set $SAMPLE
		$bbduk -Xmx1g in1=${3}_R1_demux2.fastq.gz in2=${3}_R2_demux2.fastq.gz out1=./${Run}/1_trimmed/${3}_trim_R1.fastq.gz out2=./${Run}/1_trimmed/${3}_trim_R2.fastq.gz ref=./Scripts/${Run}/barcodes/${i}/${3}_adapters.fa ktrim=r k=20 mink=5 hdist=1 overwrite=TRUE stats=./${Run}/1_trimmed/${3}_trimReport tpe tbo
	done < <(grep "	${i}	" ./Scripts/${Run}/${Run}.txt)
done


##Index reference genome
java -Xmx4g -jar CreateSequenceDictionary.jar R=panTro6.fa OUTPUT=panTro6.dict
bwa-0.7 index panTro6.fa


##Map to reference genome
#Follows a modified version of GATK vest pracice guidelines as of 26.09.2018
#Convert fastq files to uBAM files. Then to map to chimp reference genome, then clean up using MergeAlignedBam
#Requires indexed reference genome in fasta format and a list of all samples per sequencing Run
Run="Sequencing Run Name"
REF="pT6"
ref="panTro6.fa"
while read SAMPLE; do
	set $SAMPLE
	#Create read group name
	id=`zcat ${3}_trim_R1.fastq.gz | head -n1 | cut -d'@' -f2 | cut -d':' -f1`
	pu=`zcat ${3}_trim_R1.fastq.gz | head -n1 | cut -d'@' -f2 | cut -d':' -f3`
	RG="'@RG\tID:${id}\tPL:Illumina\tPU:${pu}\tLB:${8}-${10}\tSM:${4}'"
	rg="@RG\tID:${id}\tPL:Illumina\tPU:${pu}\tLB:${8}-${10}\tSM:${4}"
	echo ${RG}
	#Convert two fastq files to unaligned bam file (ubam) and then back to fq (Can delete these intermediate files later)
	java -Xmx4g -jar FastqToSam.jar FASTQ=${3}_trim_R1.fastq.gz FASTQ2=${3}_trim_R2.fastq.gz OUTPUT=${3}_trim.bam READ_GROUP_NAME=${RG} SAMPLE_NAME=${4} 
	java -Xmx4g -jar SamToFastq.jar I=${3}_trim.bam FASTQ=${3}_trim.fq INTERLEAVE=true NON_PF=false
	#Map reads
	bwa-0.7 mem -M -t 5 -p ${ref} ${3}_trim.fq > ../2_mapped/${3}_${REF}_mapped.sam
	java -Xmx4g -jar MergeBamAlignment.jar ALIGNED_BAM=../2_mapped/${3}_${REF}_mapped.sam UNMAPPED_BAM=${3}_trim.bam OUTPUT=../2_mapped/${3}_${REF}_mapped.bam CREATE_INDEX=true R=${ref}  PAIRED_RUN=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=1 ATTRIBUTES_TO_RETAIN=XS
	fi
done < ./Scripts/${Run}/${Run}.txt



##Merge libraries sequenced across different runs
#Libraries were either sequenced across two or four lanes
#Requires a list of libraries to merge (list excludes libraries identified as potentially contaminated)
while read MERGEME; do
	set $MERGEME
	test=`echo ${3} |tr -cd : | wc -c`
	#Merge bams of each library into single bam
	if [[ "$test" == 2 ]]
     	then
		lib1=`echo ${3} | cut -d"," -f1 | sed 's/:/\/2_mapped\//g'`
		lib2=`echo ${3} | cut -d"," -f2 | sed 's/:/\/2_mapped\//g'`
		echo "2 Libraries to Merge: ${lib1} ${lib2}"
		java -Xmx4g -jar MergeSamFiles.jar I=${lib1}_${REF}_mapped.bam I=${lib2}_${REF}_mapped.bam O=MergedLibs/${1}_${REF}_mapped.bam AS=TRUE SO=coordinate
	else	
		lib1=`echo ${3} | cut -d"," -f1 | sed 's/:/\/2_mapped\//g'`
		lib2=`echo ${3} | cut -d"," -f2 | sed 's/:/\/2_mapped\//g'`
		lib3=`echo ${3} | cut -d"," -f3 | sed 's/:/\/2_mapped\//g'`
		lib4=`echo ${3} | cut -d"," -f4 | sed 's/:/\/2_mapped\//g'`
		echo "4 Libraries to Merge ${lib1} ${lib2} ${lib3} ${lib4}"
		java -Xmx4g -jar MergeSamFiles.jar I=${lib1}_${REF}_mapped.bam I=${lib2}_${REF}_mapped.bam I=${lib3}_${REF}_mapped.bam I=${lib4}_${REF}_mapped.bam O=MergedLibs/${1}_${REF}_mapped.bam AS=TRUE SO=coordinate
	fi
done < ./Scripts/Merged/MergeLibs.txt


##Filter reads and mark and remove duplicates
for i in *mapped.bam; do
	echo ${i}
	out=`echo ${i} | cut -d"." -f1`
	lib=`echo ${i} | cut -d"_" -f1`
	echo "${lib} ${out}"
	MarkDuplicates.jar I=${i} O=${out}_DeDup.bam M=${lib}_MarkDup_Report.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE 
	bam-mangle -e 'MAPQ>=30 && LENGTH>35' -o ${out}_DeDup_Filtered.bam ${out}_DeDup.bam
done


##Merge individuals sequenced across different libraries
#Requires a list of libraries to merge (list excludes libraries and individuals identified as potentially contaminated)
while read MERGEME; do
	set $MERGEME
	echo ${1}
	test=${3}
	echo ${test}
	#Merge bams of each library into single bam per individual
	if [[ "$test" == 1 ]]
	then
		lib1="${2}_${REF}_mapped_DeDup_Filtered.bam"
		echo "copying ${2}"
		cp ${lib1} ../MergedInds/${1}_${REF}.bam
	fi
	if [[ "$test" == 2 ]]
     	then
		lib1=`echo ${2} | cut -d"," -f1`
		lib2=`echo ${2} | cut -d"," -f2`
		echo "2 Libraries to Merge: ${lib1} ${lib2}"
		samtools merge -f ../MergedInds/${1}_${REF}.bam ${lib1}_${REF}_mapped_DeDup_Filtered.bam ${lib2}_${REF}_mapped_DeDup_Filtered.bam
	fi	
	if [[ "$test" == 3 ]]
    	then
		lib1=`echo ${2} | cut -d"," -f1`
		lib2=`echo ${2} | cut -d"," -f2`
		lib3=`echo ${2} | cut -d"," -f3`
		echo "3 Libraries to Merge: ${lib1} ${lib2} ${lib3}"
		samtools merge -f ../MergedInds/${1}_${REF}.bam ${lib1}_${REF}_mapped_DeDup_Filtered.bam ${lib2}_${REF}_mapped_DeDup_Filtered.bam ${lib3}_${REF}_mapped_DeDup_Filtered.bam
	fi
	if [[ "$test" == 4 ]]
      	then
		lib1=`echo ${2} | cut -d"," -f1`
		lib2=`echo ${2} | cut -d"," -f2`
		lib3=`echo ${2} | cut -d"," -f3`
		lib4=`echo ${2} | cut -d"," -f4`
		echo "4 Libraries to Merge: ${lib1} ${lib2} ${lib3} ${lib4}"
		samtools merge -f ../MergedInds/${1}_${REF}.bam ${lib1}_${REF}_mapped_DeDup_Filtered.bam ${lib2}_${REF}_mapped_DeDup_Filtered.bam ${lib3}_${REF}_mapped_DeDup_Filtered.bam ${lib4}_${REF}_mapped_DeDup_Filtered.bam
	fi
	if [[ "$test" == 5 ]]
      	then
		lib1=`echo ${2} | cut -d"," -f1`
		lib2=`echo ${2} | cut -d"," -f2`
		lib3=`echo ${2} | cut -d"," -f3`
		lib4=`echo ${2} | cut -d"," -f4`
		lib5=`echo ${2} | cut -d"," -f5`
		echo "5 Libraries to Merge: ${lib1} ${lib2} ${lib3} ${lib4} ${lib5}"
		samtools merge -f ../MergedInds/${1}_${REF}.bam ${lib1}_${REF}_mapped_DeDup_Filtered.bam ${lib2}_${REF}_mapped_DeDup_Filtered.bam ${lib3}_${REF}_mapped_DeDup_Filtered.bam ${lib4}_${REF}_mapped_DeDup_Filtered.bam ${lib5}_${REF}_mapped_DeDup_Filtered.bam
	fi
	if [[ "$test" == 6 ]]
      	then
		lib1=`echo ${2} | cut -d"," -f1`
		lib2=`echo ${2} | cut -d"," -f2`
		lib3=`echo ${2} | cut -d"," -f3`
		lib4=`echo ${2} | cut -d"," -f4`
		lib5=`echo ${2} | cut -d"," -f5`
		lib6=`echo ${2} | cut -d"," -f6`
		echo "6 Libraries to Merge: ${lib1} ${lib2} ${lib3} ${lib4} ${lib5} ${lib6}"
		samtools merge -f ../MergedInds/${1}_${REF}.bam ${lib1}_${REF}_mapped_DeDup_Filtered.bam ${lib2}_${REF}_mapped_DeDup_Filtered.bam ${lib3}_${REF}_mapped_DeDup_Filtered.bam ${lib4}_${REF}_mapped_DeDup_Filtered.bam ${lib5}_${REF}_mapped_DeDup_Filtered.bam ${lib6}_${REF}_mapped_DeDup_Filtered.bam
	fi
	if [[ "$test" == 7 ]]
     	then
		lib1=`echo ${2} | cut -d"," -f1`
		lib2=`echo ${2} | cut -d"," -f2`
		lib3=`echo ${2} | cut -d"," -f3`
		lib4=`echo ${2} | cut -d"," -f4`
		lib5=`echo ${2} | cut -d"," -f5`
		lib6=`echo ${2} | cut -d"," -f6`
		lib7=`echo ${2} | cut -d"," -f7`
		echo "7 Libraries to Merge: ${lib1} ${lib2} ${lib3} ${lib4} ${lib5} ${lib7}"
		samtools merge -f ../MergedInds/${1}_${REF}.bam ${lib1}_${REF}_mapped_DeDup_Filtered.bam ${lib2}_${REF}_mapped_DeDup_Filtered.bam ${lib3}_${REF}_mapped_DeDup_Filtered.bam ${lib4}_${REF}_mapped_DeDup_Filtered.bam ${lib5}_${REF}_mapped_DeDup_Filtered.bam ${lib6}_${REF}_mapped_DeDup_Filtered.bam ${lib7}_${REF}_mapped_DeDup_Filtered.bam
	fi
done </home/lauren_white/Scripts/Process/Enriched/Merged/MergeInds.txt


##ID SNPs from high quality samples and estimate allele frequencies
#Requries a bamlist of the highquality individuals
ANGSD=path to angsd
minind=225 #90% of 10x dataset
maf=0.05 
MAF="05"
mD=675 #average of 3 reads per individual
for i in `cat Autosomes.txt`
do 
	echo ${i}
	$ANGSD -P 1 -b GeoRef10xInds_DeContam.bamlist -r ${i} -out GeoRef10xInds_DeContam_minInd${minind}_maf${MAF}_minDepth${mD}_${i} -gl 2 -doMajorMinor 1 -doMaf 1 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -rmTriallelic 1e-6 -SNP_pval 1e-6 -minInd ${minind} -minmaf ${maf} -setMinDepth ${mD} -doHWE 1 -minHWEpval 0.001 -doCounts 1
	zcat GeoRef10xInds_DeContam_minInd${minind}_maf${MAF}_minDepth${mD}_${i}.mafs.gz | cut -f1-4 | grep -v "chromo" > GeoRef10xInds_DeContam_minInd${minind}_maf${MAF}_minDepth${mD}_${i}_sites.txt
	$ANGSD sites index GeoRef10xInds_DeContam_minInd${minind}_maf${MAF}_minDepth${mD}_${i}_sites.txt
done
#Extract MAFs
mafs=`ls GeoRef10xInds*.mafs.gz`
zcat ${mafs} | cut -f5 | sed 1d | grep -v "knownEM" > GeoRef10xInds_DeContam_minInd${minind}_maf${MAF}_minDepth${mD}_Autosomes.freq


##Create genotype liklihood files IN BEAGLE BINARY FORMAT for sites ID's from GeoREf10x dataset for all individual -> For ngsRelate
#requries a list of bamfiles of all individuals and a list of SNPs (sites) ID'd from the high quality dataset above
ANGSD=path to angsd
minind=225 #90% of 10x dataset
maf=0.05 
MAF="05"
mD=675 #average of 3 reads per individual
for i in `cat Autosomes.txt`
do
	echo ${i}
	$ANGSD -P 1 -b AllInds_DeContam.bamlist -r ${i} -out ./GLF_3/AllInds_DeContam_GeoRef10xSites_minInd${minind}_maf${MAF}_minDepth${mD}_${i} -GL 2 -doGlf 3 -doMajorMinor 3 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -sites GeoRef10xInds_DeContam_minInd${minind}_maf${MAF}_minDepth${mD}_${i}_sites.txt 
done


##Use Genotype liklihoods to estimate relatedness 
ANGSD=path to angsd
ngsr=path to ngsRelate
minind=225 #90% of 10x dataset
maf=0.05 
MAF="05"
mD=675 #average of 3 reads per individual
#Put SNPs from the autosomes into one file
beag=`ls AllInds_DeContam_GeoRef10xSites*.glf.gz`
zcat ${beag} > ./GLF_3/AllInds_DeContam_GeoRef10xSites_minInd${minind}_maf${MAF}_minDepth${mD}_Autosomes.glf
gzip ./GLF_3/AllInds_DeContam_GeoRef10xSites_minInd${minind}_maf${MAF}_minDepth${mD}_Autosomes.glf
$ngsr -p 8  -g ./GLF_3/AllInds_DeContam_GeoRef10xSites_minInd${minind}_maf${MAF}_minDepth${mD}_Autosomes.glf.gz -n 482 -f GeoRef10xInds_DeContam_minInd${minind}_maf${MAF}_minDepth${mD}_Autosomes.freq -z AllInds_DeContam_SampleNames.txt -O AllInds_DeContam_GeoRef10xSites_minInd${minind}_maf${MAF}_minDepth${mD}_GL_Autosomes_Relate.txt


