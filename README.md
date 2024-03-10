# STAR+WASP reduces reference bias in the allele-specific mapping of RNA-seq reads
This repository contains workflows, analytics files and code used in the generation of the manuscript: *STAR+WASP reduces reference bias in the allele-specific mapping of RNA-seq reads*  
<!--- &nbsp; --->
##### File descriptions:
-------------------------

| **SN** |**Directory** | **File**   | **Description** |
|----------------|------------|------------|------------|
|1|[read_alignment](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/read_alignment)|[base_script_STAR_runs.sh](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/read_alignment/base_script_STAR_runs.sh)|This file examplifies how read alignments with STAR were conducted|
|2|[read_alignment](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/read_alignment)|[base_script_STAR_WASP_runs.sh](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/read_alignment/base_script_STAR_WASP_runs.sh)|This file examplifies how read alignments with STAR+WASP were conducted|
|3|[read_alignment](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/read_alignment)|[base_script_WASP_runs.sh](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/read_alignment/base_script_WASP_runs.sh)|This file examplifies how read alignments with WASP were conducted|
|4|[read_alignment](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/read_alignment)|[Prep_STAR_Runs.py](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/read_alignment/Prep_STAR_Runs.py)|STAR runs - prep file|
|5|[read_alignment](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/read_alignment)|[Prep_STAR_WASP_Runs.py](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/read_alignment/Prep_STAR_WASP_Runs.py)|STAR+WASP runs - prep file|
|6|[read_alignment](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/read_alignment)|[Prep_WASP_Runs.py](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/read_alignment/Prep_WASP_Runs.py)|WASP runs - prep file|
|7|[downstream_analytics](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/downstream_analytics)|[benchmark_results_all_runs.txt](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/downstream_analytics/benchmark_results_all_runs.txt)|This file contains extracted benchmak results for all runs to encapsulate the wall clock per alignment, system and user time, maximum memory used per run, e.t.c. for runs conducted on a dedicated server|
|8|[downstream_analytics](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/downstream_analytics)|[benchmark_results_shared_computing_environment.txt](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/downstream_analytics/benchmark_results_shared_computing_environment.txt)|This file contains extracted benchmak results for all runs to encapsulate the wall clock per alignment, start and end time, maximum memory used per run, e.t.c. for runs conducted in a shared computing environment|
|9|[downstream_analytics](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/downstream_analytics)|[comp_runs.sh](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/downstream_analytics/comp_runs.sh)|Bash file used to extract data from STAR+WASP and WASP alignments for comparative analyses|
|10|[downstream_analytics](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/downstream_analytics)|[downstream_analytics.Rmd](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/downstream_analytics/downstream_analytics.Rmd)|RMD file that shows data analysis and visualization used to generate manuscript figures|
|11|[downstream_analytics](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/downstream_analytics)|[final_log_results_all_runs.txt](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/downstream_analytics/final_log_results_all_runs.txt)|This ......|
|12|[downstream_analytics](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/tree/main/downstream_analytics)|[final_log_results_extract_all_runs.py](https://github.com/rasiimwe/STAR-WASP-reduces-reference-bias-in-the-allele-specific-mapping-of-RNA-seq-reads/blob/main/downstream_analytics/final_log_results_extract_all_runs.py)|This ......|


<!--- --->

&nbsp;
&nbsp;

-------------------------
**Please cite (subject to change):**

Rebecca Asiimwe, Dobin Alexander. STAR+WASP reduces reference bias in the allele-specific mapping of RNA-seq reads. bioRxiv 2024:2024.01.21.576391. (https://www.biorxiv.org/content/10.1101/2024.01.21.576391v2)

doi: https://doi.org/10.1101/2024.01.21.576391
