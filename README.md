# pileup2vaf
detect snv and small indel using pysam pileup method

# Usage
`python3 pileup2vaf.py -bam $bam -bed $bed -fa $fa -od $PWD -n "rawbam"`

```
usage: Detect snv and small indel using pysam pileup method
       [-h] [-bam BAM] [-bed BED] [-d DEPTH] [-n NAME] [-mq MAPQ] [-bq BASEQ]
       [-fa FASTA] [-od OUTDIR]

optional arguments:
  -h, --help  show this help message and exit
  -bam BAM    bam file
  -bed BED    bed file
  -d DEPTH    max depth
  -n NAME     samplel name
  -mq MAPQ    map quality
  -bq BASEQ   base quality
  -fa FASTA   fasta file
  -od OUTDIR  output dir
```
`-bam`: input bam file

`-bed`: BED file

`-fa` : reference fasta file

`-d`  : max depth limit, same as samtools mpileup. [default: 8000]

`-n`  : sample name

`-mq` : mapping quailty cutoff. [default: 0]

`-bq` : base quality. [default: 0]

`-od` : output dir

# Test
`less test/cmd.sh`

```
python3 ../pileup2vaf.py -bam $line -bed $bed -fa $fa -od $PWD -n $name
```


# Outfiles
The outfile is `<sample_name>.variants.xls`

## Outfile Format
```
Chr     Pos             Ref         Alt          AltNum      Depth      AltAlleleFrequency
chr12   52380364        T           ref-homo     0           24         0.0     
chr12   52380365        T           C            1           24         0.042
chr12   52380366        A           ref-homo     0           24         0.0
chr12   52380623        TGATGA      T            90          582        0.155
chrX    48649699        C           CTACT        30          155        0.194
```


# Dependency
* Python3
* Pysam
