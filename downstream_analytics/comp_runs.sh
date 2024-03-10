#Extracting data for comparative analyses from STAR_WASP and WASP runs (example extract from 32-thread runs)

#Sorting and creating indexes for bam files
samtools sort A_sorted.keep.bam -o A_sorted.keep.bam
samtools index A_sorted.keep.bam A_sorted.keep.bai

samtools sort A_sorted.to.remap.bam -o A_sorted.to.remap.bam
samtools  index  A_sorted.to.remap.bam A_sorted.to.remap.bai 


#Extracting vW tagged reads from STAR+WASP runs (reminder: no tags in STAR or WASP runs)
#Extracting overall number of reads in bamfile
samtools view -c /home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/sample/32threads/Aligned.sortedByCoord.out.bam

samtools view /home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/sample/32threads/Aligned.sortedByCoord.out.bam | grep vW:i> STAR_vW_Tagged_Reads 
wc -l STAR_vW_Tagged_Reads 
echo "Overall number of reads in .bam file:" >> comp_res.txt 
samtools view -c  /home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/sample/32threads/Aligned.sortedByCoord.out.bam >> comp_res.txt 
echo " " >> comp_res.txt
echo "STAR_vW_Tagged_Reads:" >> comp_res.txt 
wc -l STAR_vW_Tagged_Reads >> comp_res.txt 
echo " " >> comp_res.txt


#Extracting unique read-tag pairings for mapping purposes downstream
samtools view /home/asiimwe/projects/run_env/alpha_star_wasp_benchmarking/STAR_WASP/STAR_WASP_Runs/sample/32threads/Aligned.sortedByCoord.out.bam | grep vW:i | awk '{print $1 "\t" $(NF)}' | sort -u>  STAR_vW_Tagged_Reads_Unique_Subset
awk -F, 'a[$1]++{print $1}' STAR_vW_Tagged_Reads_Unique_Subset #rechecking for duplicates
wc -l STAR_vW_Tagged_Reads
wc -l STAR_vW_Tagged_Reads_Unique_Subset
echo "STAR_vW_Tagged_Reads_Unique_Subset:" >> comp_res.txt 
wc -l STAR_vW_Tagged_Reads_Unique_Subset >> comp_res.txt 
echo " " >> comp_res.txt

#Mapping WASP reads to STAR vW tagged reads by intersecting the two files to associate WASP reads with their respective tags
#1. First, we will sort and extract all unique WASP reads
samtools view -c A_sorted.bam
echo "WASP - Overall number of reads in .bam file:" >> comp_res.txt 
samtools view -c A_sorted.bam >> comp_res.txt 
echo " " >> comp_res.txt

samtools view A_sorted.bam | sort -u -k1,1 > WASP_Reads_Sorted_Unique 
#Checking for duplicates:
awk -F, 'a[$1]++{print $1}' WASP_Reads_Sorted_Unique
#samtools view -c A_sorted.bam 
echo "WASP_Reads_Sorted_Unique:" >> comp_res.txt
wc -l WASP_Reads_Sorted_Unique  >> comp_res.txt 
echo " " >> comp_res.txt

#2. Next we'll associate WASP-run reads with the vW tags from the STAR+WASP run by joining the "WASP_Reads_Sorted_Unique" file with "STAR_vW_Tagged_Reads_Unique_Subset (this file contains only 2 columns, read and associated tag - join will be on read id)" 
#If the files are not sorted, run the following before the join:
sort -u -k1,1 STAR_vW_Tagged_Reads_Unique_Subset > Sorted_STAR_vW_Tagged_Reads_Unique_Subset
mv Sorted_STAR_vW_Tagged_Reads_Unique_Subset STAR_vW_Tagged_Reads_Unique_Subset
join -1 1 -2 1 WASP_Reads_Sorted_Unique STAR_vW_Tagged_Reads_Unique_Subset > STAR_WASP_vW_Tagged_Reads_Unique
wc -l STAR_WASP_vW_Tagged_Reads_Unique
awk -F, 'a[$1]++{print $1}' STAR_WASP_vW_Tagged_Reads_Unique
echo "STAR_WASP_vW_Tagged_Reads_Unique:" >> comp_res.txt
wc -l STAR_WASP_vW_Tagged_Reads_Unique  >> comp_res.txt 
echo " " >> comp_res.txt

#Extracting only read-tag pairs
awk '{print $1 "\t" $(NF)}' STAR_WASP_vW_Tagged_Reads_Unique >  STAR_WASP_vW_Tagged_Reads_Unique_Subset


#Checking WASP reads that need remapping:
echo "Overall number of reads to remap:" >> comp_res.txt
samtools view -c A_sorted.to.remap.bam >> comp_res.txt 
echo " " >> comp_res.txt

samtools view A_sorted.to.remap.bam | awk '{print $1}' | sort -u> reads_to_remap_unique
echo "Unique reads to remap:" >> comp_res.txt
wc -l reads_to_remap_unique >> comp_res.txt #These are WASP reads that need remapping
echo " " >> comp_res.txt

#This is overall file with remapped reads, we want to extract vW_tagged reads so we shall join these files as well to provide that mapping
join -1 1 -2 1 reads_to_remap_unique STAR_vW_Tagged_Reads_Unique_Subset > reads_to_remap_unique_tagged

#WASP reads to keep
echo "Overall number of reads to keep:" >> comp_res.txt
samtools view -c A_sorted.keep.bam >> comp_res.txt 
echo " " >> comp_res.txt

samtools view A_sorted.keep.bam |awk '{print $1}' | sort -u> reads_to_keep_unique
echo "Unique WASP reads to keep:" >> comp_res.txt
wc -l reads_to_keep_unique >> comp_res.txt
echo " " >> comp_res.txt
join -1 1 -2 1 reads_to_keep_unique STAR_vW_Tagged_Reads_Unique_Subset > reads_to_keep_unique_tagged


#Comparative matrix files: STAR_vW_Tagged_Reads_Unique_Subset(unique reads with tags), reads_to_remap_unique, reads_to_keep_unique
#Checking duplicates in all input files"
awk -F, 'a[$1]++{print $1}' STAR_WASP_vW_Tagged_Reads_Unique_Subset
awk -F, 'a[$1]++{print $1}' reads_to_remap_unique_tagged
awk -F, 'a[$1]++{print $1}' reads_to_keep_unique_tagged
#All files have no duplicate reads - Consider including branching logic for instances where we may have or need to consider duplications


awk 'FNR==NR{a[$1]=$1;next}{print $0,a[$1]?a[$1]:"NA"}' reads_to_remap_unique_tagged STAR_WASP_vW_Tagged_Reads_Unique_Subset > pass1_to_remap.txt # we can either tag reads that do not intersect as NA or 0: {print $0,a[$1]?a[$1]:0}
#We also expect wc -l for both pass1_to_remap and pass2 below == wc -l STAR_WASP_vW_Tagged_Reads_Unique_Subset. We are takig all reads in STAR_WASP_vW_Tagged_Reads_Unique_Subset and seeing which of those were flagged for remapping or keeping. If true, fill = read id, if false, fill = NA

#Next we take the output from pass1 and do the same mapping for reads that were kept - column concatenated. Also if true, fill = read id, if false, fill = NA
awk 'FNR==NR{a[$1]=$1;next}{print $0,a[$1]?a[$1]:"NA"}' reads_to_keep_unique_tagged pass1_to_remap.txt > pass2_to_keep.txt

#Extracting file with read, chr, pos, vA, vG and vW_tags for reference bias analysis
awk '{print $1, $3, $4, $22, $23, $24}' STAR_vW_Tagged_Reads | sort -u -k1,1  > STAR_vW_Tagged_Reads_vA_vG

