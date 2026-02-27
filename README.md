# Differential expression and functional analysis of yeast transcriptomics data


## **1 | Introduction**


## **2 | Methods**

## **2.1 Data Description**
```
prefetch
fasterq-dump

datasets download genome accession GCF_000146045.2 --include genome,rna,gtf --filename yeast_ref.zip
```
Reference genome, transcriptome, and annotation file renamed to ref.fna, rna.fna, and genomic.gtf, respectively

Creating decoy file
```
grep "^>" ref/ref.fna | cut -d " " -f 1 > decoys.txt
sed -i 's/>//g' decoys.txt
```

Combined genome + transcriptome
```
cat ref/rna.fna ref/ref.fna > genome_transcripts.fa

salmon index -t genome_transcripts.fa -d decoys.txt -p 6 -i transcripts_index
```

## **3 | Results**


## **4 | Discussion**


## References
