# popscle
popscle is a suite of population scale analysis tools for single-cell genomics data including implementation of Demuxlet / Freemuxlet methods and auxilary tools

### Introduction

`demuxlet` and `freemuxlet` are two software tools to deconvolute sample identity and identify multiplets when multiple samples are pooled by barcoded single cell sequencing. If external genotyping data for each sample is available, `demuxlet` would be recommended. On the other hand, if external genotyping data is not available, the genotyping-free version demuxlet, `freemuxlet`, would be recommended. 

To reduce computation time, we need to run `dsc-pileup` before running `demuxlet` and `freemuxlet`. `dsc-pileup` is a software tool to pileup reads and corresponding base quality for  each overlapping SNPs and each barcode. By using pileup files, it would allow us to run demuxlet/freemuxlet pretty fast multiple times without going over the BAM file again. 

`dsc-pileup` requires the following input files:

1. a SAM/BAM/CRAM file produced by the standard 10x sequencing platform, or any other barcoded single cell RNA-seq (with proper `--tag-UMI` and `--tag-group`) options.
2. A VCF/BCF files containing (AC) and (AN) from referenced population (e.g. 1000g).

`demuxlet` require the following input files:

1. Pileup files (CEL,VAR and PLP) produced by `dsc-pileup`.
2. a VCF/BCF file containing the genotype (GT), posterior probability (GP), or genotype likelihood (GL) to assign each barcode to a specific sample (or a pair of samples) in the VCF file.

Alternatively, `dsc-pileup` could also directly take SAM file without running `dsc-pileup`. In this case, `dsc-pileup` would require the following files:

1. a SAM/BAM/CRAM file produced by the standard 10x sequencing platform, or any other barcoded single cell RNA-seq (with proper `--tag-UMI` and `--tag-group`) options.
2. a VCF/BCF file containing the genotype (GT), posterior probability (GP), or genotype likelihood (GL) to assign each barcode to a specific sample (or a pair of samples) in the VCF file.

`freemuxlet` require the following input: 

1. Pileup files (CEL, PLP and VAR) from dsc-pileup 
2. Number of samples

### Tips for running
* If external reference sequence vcf file is available, **_demuxlet_** is recommended
* Default setting alpha as 0.5, which assumes the expected proportion of 50% genetic mixture from two individuals, to get better estimates of doublets.
* Set `--group-list` to a list of barcodes (i.e. barcodes.tsv from 10X) in `dsc-pileup` to speed things up and only get demultiplexing for cells called by other methods.
* To reproduce the results presented in Figure 2 of the demuxlet paper, please go to: https://github.com/yelabucsf/demuxlet_paper_code/tree/master/fig2 to download the vcf and the outputs of demuxlet.
* Check tutorial README.md for more detailed tutorial with example data

### Installing demuxlet/freemuxlet


$ mkdir build

$ cd build

$ cmake ..

<pre>
In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/  -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.a ..

For bzip2:
  - $ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - $ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
</pre>

$ make


### Using demuxlet and freemuxlet
All softwares use a self-documentation utility. You can run each utility with -man or -help option to see the command line usages. Also, we offer some general practice with an example in tutorial (data is available here: https://ucsf.box.com/s/vg1bycvsjgyg63gkqsputprq5rxzjl6k).

#### demuxlet
<pre>
$(POPSCLE_HOME)/bin/popscle dsc-pileup --sam /data/$bam --vcf /data/$ref_vcf --out /data/$pileup
$(POPSCLE_HOME)/bin/popscle demuxlet --plp /data/$pileup --vcf /data/$external_vcf --field $(GT or GP or PL) --out /data/$filename
</pre>

Or, demuxlet could directly take SAM file as input:
<pre>
$(POPSCLE_HOME)/bin/popscle demuxlet --sam /data/$sam --vcf /data/$external_vcf --field $(GT or GP or PL) --out /data/$filename
</pre>

#### freemuxlet
<pre>
$(POPSCLE_HOME)/bin/popscle dsc-pileup --sam /data/$bam --vcf /data/$ref_vcf --out /data/$pileup
$(POPSCLE_HOME)/bin/popscle freemuxlet --plp /data/$pileup --out /data/$filename --nsample $n
</pre>

The detailed usage is also pasted below.

#### dsc-pileup
<pre>
Options for input SAM/BAM/CRAM 
   --sam [STR: ] : Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed
   --tag-group [STR: CB] : Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB
   --tag-UMI [STR: UB] : Tag representing UMIs. For 10x genomiucs, use UB

Options for input VCF/BCF
   --vcf [STR: ] : Input VCF/BCF file, containing the AC and AN field
   --sm [V_STR: ] : List of sample IDs to compare to (default: use all)
   --sm-list [STR: ] : File containing the list of sample IDs to compare

Output Options
   --out [STR: ]  : Output file prefix
   --sam-verbose [INT: 1000000] : Verbose message frequency for SAM/BAM/CRAM
   --vcf-verbose [INT: 10000] : Verbose message frequency for VCF/BCF
   --skip-umi [FLG: OFF] : Do not generate [prefix].umi.gz file, which stores the regions covered by each barcode/UMI pair

SNP-overlapping Read filtering Options
   --cap-BQ [INT: 40] : Maximum base quality (higher BQ will be capped)
   --min-BQ [INT: 13] : Minimum base quality to consider (lower BQ will be skipped)
   --min-MQ [INT: 20] : Minimum mapping quality to consider (lower MQ will be ignored)
   --min-TD [INT: 0] : Minimum distance to the tail (lower will be ignored)
   --excl-flag [INT: 3844] : SAM/BAM FLAGs to be excluded

Cell/droplet filtering options
   --group-list [STR: ] : List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run
   --min-total [INT: 0] : Minimum number of total reads for a droplet/cell to be considered
   --min-uniq [INT: 0] : Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered
   --min-snp [INT: 0] : Minimum number of SNPs with coverage for a droplet/cell to be considered
</pre>

#### demuxlet
<pre>
Options for input SAM/BAM/CRAM
   --sam [STR: ] : Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed
   --tag-group [STR: CB] : Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB
   --tag-UMI [STR: UB] : Tag representing UMIs. For 10x genomiucs, use UB

Options for input Pileup format
   --plp [STR: ] : Input pileup format

Options for input VCF/BCF
   --vcf [STR: ] : Input VCF/BCF file, containing the individual genotypes (GT), posterior probability (GP), or genotype likelihood (PL)
   --field [STR: GP] : FORMAT field to extract the genotype, likelihood, or posterior from
   --geno-error-offset [FLT: 0.10] : Offset of genotype error rate. [error] = [offset] + [1-offset]*[coeff]*[1-r2]
   --geno-error-coeff  [FLT: 0.00] : Slope of genotype error rate. [error] = [offset] + [1-offset]*[coeff]*[1-r2]
   --r2-info [STR: R2] : INFO field name representing R2 value. Used for representing imputation quality
   --min-mac [INT: 1] : Minimum minor allele frequency
   --min-callrate [FLT: 0.50] : Minimum call rate
   --sm [V_STR: ] : List of sample IDs to compare to (default: use all)
   --sm-list [STR: ] : File containing the list of sample IDs to compare

Output Options
   --out [STR: ] : Output file prefix
   --alpha [V_FLT: ] : Grid of alpha to search for (default is 0.1, 0.2, 0.3, 0.4, 0.5)
   --doublet-prior [FLT: 0.50] : Prior of doublet
   --sam-verbose [INT: 1000000] : Verbose message frequency for SAM/BAM/CRAM
   --vcf-verbose [INT: 10000] : Verbose message frequency for VCF/BCF

Read filtering Options
   --cap-BQ [INT: 40] : Maximum base quality (higher BQ will be capped)
   --min-BQ [INT: 13] : Minimum base quality to consider (lower BQ will be skipped)
   --min-MQ [INT: 20] : Minimum mapping quality to consider (lower MQ will be ignored)
   --min-TD [INT: 0] : Minimum distance to the tail (lower will be ignored)
   --excl-flag [INT: 3844] : SAM/BAM FLAGs to be excluded

Cell/droplet filtering options
   --group-list [STR: ] : List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run
   --min-total [INT: 0] : Minimum number of total reads for a droplet/cell to be considered
   --min-uniq [INT: 0] : Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered
   --min-snp [INT: 0] : Minimum number of SNPs with coverage for a droplet/cell to be considered
</pre>

#### freemuxlet
<pre>
Options for input pileup
--plp [STR: ] : Prefix of input files generated by dsc-pileup
--init-cluster [STR: ] : Input file containing the initial cluster information

Output Options
--out [STR: ] : Output file prefix
--nsample [INT: 0] : Number of samples multiplexed together
--aux-files [FLG: OFF] : Turn on writing auxiliary output files
--verbose [INT: 100] : Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages

Options for statistical inference
--doublet-prior [FLT: 0.50] : Prior of doublet
--bf-thres [FLT: 5.41] : Bayes Factor Threshold used in the initial clustering
--frac-init-clust [FLT: 0.50] : Fraction of droplets to be clustered in the very first round of initial clustering procedure
--iter-init [INT: 10] : Iteration for initial cluster assignment (set to zero to skip the iterations)
--keep-init-missing [FLG: OFF] : Keep missing cluster assignment as missing in the initial iteration

Read filtering Options
--cap-BQ [INT: 40] : Maximum base quality (higher BQ will be capped)
--min-BQ [INT: 13] : Minimum base quality to consider (lower BQ will be skipped)

Cell/droplet filtering options
--group-list [STR: ] : List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run
--min-total [INT: 0] : Minimum number of total reads for a droplet/cell to be considered
--min-uniq [INT: 0] : Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered
--min-snp [INT: 0] : Minimum number of SNPs with coverage for a droplet/cell to be considered
</pre>

### Interpretation of output files

#### dsc-pileup
**_dsc-pileup_** generates multiple output file, such as  `[prefix].cel`, `[prefix].var`, `[prefix].plp` and `[prefix].umi`, and these three files would be the input for **_freemuxlet_** and **_demuxlet_**.

* The `[prefix].cel` file contains the relation between numerated barcode ID and barcode. Also, it contains the number of SNP and number of UMI for each barcoded droplet.  
* The `[prefix].plp` file contains the overlapping SNP and the corresponding read and base quality for each barcode ID.
* The `[prefix].var` file contains the position, reference allele and allele frequency for each SNP. 
* The `[prefix].umi` file contains the position covered by each umi.


#### demuxlet/freemuxlet
Both **_demuxlet_** and **_freemuxlet_** generate output file contains the best guess of the sample identity, with detailed statistics to reach to the best guess. It is called  `[prefix].best` for **_demuxlet_** and  `[prefix].clust1.samples.gz` for **_freemuxlet_** 

Both `[prefix].best` and `[prefix].clust1.samples.gz` file contains the following 19 columns, but `[prefix].clust1.samples.gz` contains one additional INT_ID column showing the numerated BARCODE ID. 

 1. BARCODE - Cell barcode for the cell that is being assigned in this row
 2. RD.TOTL - The total number of reads overlapping with variant sites for each droplet.
 3. RD.PASS - The total number of reads that passed the quality threshold, such as mapping quality, base quality. 
 4. RD.UNIQ - The total number of UMIs that passed the quality threshold. If a UMI is observed in a single variant multiple times, it won't be counted more. If a UMI is observed across multiple variants, it will be counted as different.
 5. N.SNP   - The total number of variants overlapping with any read in the droplet.
 6. BEST    - The best assignment for sample ID.
    * For singlets, SNG-<sample ID>
    * For doublets, DBL-<sample ID1>-<sampleID2>-<mixture rate>
    * For ambiguous droplets, , AMB-<best-singlet-sampleID>-<next-best-singlet-sampleID>-<doublet ID1/ID2>)
 7. SNG.1ST - The best singlet assignment for sample ID
 8. SNG.LLK1 - The log(likelihood that the ID from SNG.1ST is the correct assignment)       
 9. SNG.2ND - The next best singlet assignment for sample ID
 10. SNG.LLK2 - The log(likelihood that the ID from SNG.2ND is the correct assignment)        
 11. SNG.LLK0 - The log-likelihood from allele frequencies only      
 12. DBL.1ST - The sample ID that is most likely included if the assignment is a doublet
 13. DBL.2ND - The sample ID that is next most likely included ifthe assignment is a doublet
 14. ALPHA   - % Mixture Proportion
 15. LLK12   - The log(likelihood that the ID is a doublet)
 16. LLK1    - The log(likelihood that the ID from DBL.1ST is the correct singlet assignment)
 17. LLK2    - The log(likelihood that the ID from DBL.2ND is the correct singlet assignment)
 18. LLK10   - The log(likelihood that the ID from DBL.1ST is one of the doublet, and the other doublet identity is calculated from allele frequencies only)   
 19. LLK20   - The log(likelihood that the ID from DBL.2ND is one of the doublet, and the other doublet identity is calculated from allele frequencies only)   
 20. LLK00   - The log(likelihood that the droplet is doublet, but both identities are calculated from allele frequencies only)
 21. PRB.DBL - Posterior probability of the doublet assignment
 22. PRB.SNG1 - Posterior probability of the singlet assignment

#### Additional output files from freemuxlet
**_freemuxlet_** generates additional output files, such as `[prefix].clust1.vcf.gz`, `[prefix].lmix`, and optionally `[prefix].clust0.samples.gz`, `[prefix].clust0.vcf.gz` and `[prefix].ldist.gz` (with `--aux-files` argument). Each file contains the following information

* The `[prefix].clust1.vcf.gz` file is the vcf file for each sample inferred and clustered from `freemuxlet`
* The `[prefix].lmix` file contains basic statistics for each barcode
* The `[prefix].clust0.samples.gz` files contains the best sample identity assuming all droplets are singlets
* The `[prefix].clust0.vcf.gz` files is the vcf similar to `[prefix].clust1.vcf.gz` but assuming all droplets are singlets 
* The `[prefix].ldist.gz` files contains the pairwise Bayes factor for each possible pair of droplets.  
  

