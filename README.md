RTEA
====
RNA-TE fusion finder

## Dependency
1. [fastp](https://github.com/OpenGene/fastp)
2. [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
3. [samtools](http://www.htslib.org/)
4. [scallop](https://github.com/Kingsford-Group/scallop)
5. [bamtools](https://github.com/pezmaster31/bamtools)
6. [bwa](http://bio-bwa.sourceforge.net/)
7. [awk](https://www.gnu.org/software/gawk/)
7. R packages
   * magrittr
   * data.table
   * stringr
   * optparse
   * Rcpp
   * GenomicAlignments
   * BSgenome.Hsapiens.UCSC.hg19
   * BSgenome.Hsapiens.UCSC.hg38

## Usage
```bash
./rtea_pipeline <R1.fq.gz> <R2.fq.gz> <sample_name> <hisat2_ref> <threads> <out_dir> <build>
```
**hisat2_ref**: Unpacked [hisat2 index](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data) path
