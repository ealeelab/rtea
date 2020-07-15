RTEA
====
RNA-TE fusion finder

## Dependency
1. [fastp](https://github.com/OpenGene/fastp)
2. [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
3. [samtools](http://www.htslib.org/)
4. [bamtools](https://github.com/pezmaster31/bamtools)
5. [bwa](http://bio-bwa.sourceforge.net/)
6. [awk](https://www.gnu.org/software/gawk/)
7. R packages
   * magrittr
   * data.table
   * stringr
   * optparse
   * Rcpp
   * GenomicAlignments
   * BSgenome.Hsapiens.UCSC.hg19
   * BSgenome.Hsapiens.UCSC.hg38
	 * EnsDb.Hsapiens.v75
	 * EnsDb.Hsapiens.v86

## Usage
```bash
./rtea_pipeline <R1.fq.gz> <R2.fq.gz> <sample_name> <hisat2_ref> <threads> <out_dir> <build>
```
**hisat2_ref**: Unpacked [hisat2 index](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data) path


## Output
1. chr: Chromosome name
2. pos: Fusion breakpoint position on the chromosome
3. ori: Fusion direction on the chromosome (f, TE|gene; r, gene|TE)
5. class: TE class
9. seq: Proximal portion of fusion sequence
10. isPolyA: Whether it is a fusion with polyA sequence
11. posRepFamily: Repeat masked repeat family on the breakpoint position
12. posRep: Repeat maskec repeat element on the breakpoint position
16. TEfamily: TE family with highest alingment score when fusion sequence is aligned with consensus TE sequence.
17. TEscore: Alignment score of fusion sequence with the consensus TE sequence.
18. TEside: Fusion direction on the consensus TE sequence (5, TE|gene; 3, gene|TE)
19. TEbreak: Fusion breakpoint position on the consensus TE sequence.
20. depth: Number of RNA-seq reads on the breakpoint position
21. matchCnt: Number of fusion supporting RNA-seq reads.
25. polyAcnt: Number of polyA reads
28. baseQual: Median base quality of supporting reads.
29. lowMapQual: Number of supporting reads that has low mapping quality
32. mateDist: Minimum distance of mate reads
34. overhang: Distance of breakpoint from splice site
35. gap: Length of nearby intron
36. secondary: Proportion supporting reads that is from secondary alignment
38. nonspecificTE: Mean alignment score of supporting reads to consensus TE sequence
39. r1pstrand: Proportion of supporting reads that is from positive strand of chromosome
40. fusion_tx_id: Transcript ID of the fusion transcript
41. tx_support_exon: Number of read fragments spanning exonic region of the fusion transcript ID
42. tx_support_intron: Number of read gaps matching the fusion transcript ID
44. strand: Strand of fusion transcript
49. pos_type: Genomic region of breakpoint
54. polyTE: Known non-reference TE on the breakpoint position
55. hardstart: Start position of nearby reference genome where fusion sequence came from
56. hardend: Start position of nearby reference genome where fusion sequence came from
57. hardTE: Repeat masked TE subfamily of nearby reference genome where fusion sequence came from
58. hardDist: Distance from fusion breakpoint to nearby reference genome where fusion sequence came from
61. fusion_type: Type of TE fusion
62. fusion_tx_biotype: Biotype of fusion transcript
63. fusion_gene_id: Gene ID of fusion transcript
64. fusion_gene_name: Gene symbol of fusion transcript
65. Filter: Filter reason of low confidence fusion
