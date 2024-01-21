# STAR+WASP reduces reference bias in the allele-specific mapping of RNA-seq reads
This repository contains workflows, analytics files and code used in the generation of the manuscript: STAR+WASP reduces reference bias in the allele-specific mapping of RNA-seq reads  

## Read alignments: STAR, STAR+WASP and WASP
### STAR:
STAR=/usr/bin/STAR
genomeDirectory="/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/genome_directory/"
fastqFileDir="/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/FASTQ"
ulimit -n 10000

STARdir=/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR/STAR_Runs/sample_id/xthreads

cd $STARdir

STARpar="--runThreadN x --genomeDir $genomeDirectory --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD jM jI rB MC  --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c --readFilesIn $fastqFileDir/sample_id/R1.fastq.gz  $fastqFileDir/sample_id/R2.fastq.gz"

$STAR $STARpar $readFiles


### STAR+WASP 
STAR=/usr/bin/STAR
genomeDirectory="/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/genome_directory/"
fastqFileDir="/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/FASTQ"
vcfFileDir="/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/VCF"
ulimit -n 10000

STAR_WASPdir=/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/sample_id/xthreads

cd $STAR_WASPdir

STARpar="--runThreadN x --genomeDir $genomeDirectory  --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD jM jI rB MC vA vG vW --waspOutputMode SAMtag  --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c --readFilesIn $fastqFileDir/sample_id/R1.fastq.gz  $fastqFileDir/sample_id/R2.fastq.gz"

hetVcf="--varVCFfile $vcfFileDir/ssd_input_snp_dir/ssd.vcf.snv1het"


$STAR $STARpar $readFiles $hetVcf


### WASP:
WASP=/home/asiimwe/WASP
PYTHON=/home/asiimwe/miniconda3/bin/python3.9
export PATH=/usr/bin/samtools/:$PATH
STAR=/usr/bin/STAR
genomeDirectory="/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/genome_directory/"
vcfFileDir=/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/VCF
fastqFileDir=/scratch/asiimwe/STAR-WASP_FASTQs_VCFs/FASTQ
ulimit -n 10000 

WASPdir=/home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/WASP/WASP_Runs/sample_id/xthreads

cd $WASPdir

STARpar="--runThreadN x --genomeDir $genomeDirectory --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD jM jI  --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 1"

readFiles="--readFilesCommand gunzip -c --readFilesIn $fastqFileDir/sample_id/R1.fastq.gz  $fastqFileDir/sample_id/R2.fastq.gz"

hetVcf="--varVCFfile $vcfFileDir/ssd_input_snp_dir/ssd.vcf.snv1het" 

$STAR $STARpar $readFiles $hetVcf

mv $WASPdir/Aligned.sortedByCoord.out.bam $WASPdir/A_sorted.bam
samtools index $WASPdir/A_sorted.bam $WASPdir/A_sorted.bai
$PYTHON $WASP/mapping/find_intersecting_snps.py --is_paired_end --is_sorted --snp_dir $vcfFileDir/ssd_input_snp_dir/SNPdir --output_dir ./ A_sorted.bam 
$STAR $STARpar $hetVcf --readFilesCommand gunzip -c --readFilesIn  A_sorted.remap.fq1.gz A_sorted.remap.fq2.gz  
mv Aligned.sortedByCoord.out.bam Aligned.out.bam
samtools index Aligned.out.bam Aligned.out.bai
$PYTHON $WASP/mapping/filter_remapped_reads.py A_sorted.to.remap.bam  Aligned.out.bam A_sorted.keep.bam




