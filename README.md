# CFRL18_SNP_support
Scripts used to assess SNP confidence for the CFRL18 study

## SNP_support_assessment

### get_read_info_at_pos.py

Look up SNP locations in reads in .sam files and report information about the reads at those locations in a file with the following columns:

1.Contig	2.Position	3.Contig_Position	4.Reference_Reads_agree	5.Assembly_Reads_agree	6.Reference_base	7.Isolate_base	8.Assembly_base	9.Coverage	10.Consensus_score	11.Read_base_calls	12.Read_qualities	13.Read_identifiers
1. The contig in which the SNP is found.
2. The position within that contig where the SNP is found (base 1 position - i.e. starts counting at position 1).
3. The contig number and position in format Contig_Position.
4. Binary indicator of whether the reference genome and the majority base call in the reads is the same. 0 = different, 1 = same.
5. Binary indicator of whether the base called by breseq and the majority base call in the reads is the same. 0 = different, 1 = same.
6. The base in the reference genome at the position.
7. The majority base call in the reads at the position.
8. The base called by breseq at the position (if known).
9. The number of reads that align at that position (removes reads that align nearby but are clipped according to CIGAR string in .sam file.
10. Consensus score = (sum of quality scores for base calls agreeing with majority call) / (sum of all quality scores at this position among all reads).
11. Ordered list of the bases found in the reads aligned at SNP position. Order of this list corresponds to the order of the qualities and identifiers lists in columns 12 and 13.
12. Ordered list of the quality scores of base calls in reads aligned at SNP position. Quality scores are converted from phred-33 to integer form for human-readability.
13. Ordered list of read names/identifiers for reads that align to SNP position

Inputs:

**\-snp**	A file describing the locations of SNPs you want to know the support of among reads. This file should containg 1 SNP location per line and 2 tab-separated columns. The first column must be the ID of the contig in which your SNP is located. The second column must be the base-1 position (i.e. first position in the sequence is position 1)

**\-sams**	A list of paths to the .sam files you want to know the read support for base calls at your SNP locations. e.g. SNPs_directory/\*.sam

**\-ref**	Specify file with reference genome in fasta fomat (N.B. headers of contigs must match contig names in SNP file).

**\-calls** CALLS_FILE     Specify file with base calls from your assembler. File should be a tab-delimited table with a header row of SNP locations starting in the second column of format contig_position. The first column of the header line is skipped. Subsequence lines should be of the format sample_ID followed by a tab-delimitey list of bases corresponding to the position described in the header. E.g.

strain	1_142573	1_4473015
Sample_1	T	C
Sample_2	A	T

**\-outdir**	Specify the path to the directory into which you want output files written.

**\-threads**	Specify number of threads to use. Default: 1

**\-pickle**	Optional storage of reads that map to regions with snps as a pickled dict of format {contig: {1kb_bin: [list of reads]}}. This is primarily useful if you want to change parameters within the script and want to rerun the script multiple times on the same sample as it stores relevant reads and loads just those reads instead of processing the sam files every time.
