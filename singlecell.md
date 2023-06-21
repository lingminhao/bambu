# **Bambu for Single Cell Isoform Analysis**

Bambu is one of the existing bulk isoform discovery and quantification tools where a single parameter is used across samples for isoform discovery, and Bambu provides quantification that further distinguishes between full-length and non-full-length transcripts. In this work, we extend Bambu to the single-cell level, running at 12 cores whenever possible. 

## **Data** 

To use Bambu at the single-cell level, you will need the 

- **reads** the read file *(.fastq.gz) from a single-cell experiment run
- **genome** the genome file (.fasta)
- **annotations** the existing annotation file (.gtf) for the genome

## **Usage**

There are two preparatory steps before you run Bambu.

First, we use flexiplex for barcode discovery & demultiplexing. You may learn how to install flexiplex [here](https://davidsongroup.github.io/flexiplex/). Then you can run the following bash script: 

``` bash
cd flexiplex
gunzip -c ../data/GIS_cellMix_HCT116-A549-HepG2-MCF7_3primeSingleCellcDNA_Rep1_Run1.fastq.gz | ./flexiplex -p 12 -f 0
head -n 850 flexiplex_barcodes_counts.txt > my_barcode_list.txt
sort <(gunzip -c 3M-february-2018.txt.gz) <(cut -f1 my_barcode_list.txt) | uniq -d > my_filtered_barcode_list.txt
gunzip -c ../data/GIS_cellMix_HCT116-A549-HepG2-MCF7_3primeSingleCellcDNA_Rep1_Run1.fastq.gz | ./flexiplex -p 12 -k my_filtered_barcode_list.txt | gzip > ../data/new_GIS_cellMix_HCT116-A549-HepG2-MCF7_3primeSingleCellcDNA_Rep1_Run1.fastq.gz
cd ..
```
Next, we map the demultiplexed read file to the genome using minimap 2. You may learn how to install minimap2 [here](https://github.com/lh3/minimap2). 

``` bash
minimap2 -d ref.mmi data/genome
# TO DO 
```
After preparing the bam file, we can run bambu for transcript discovery and generate `readGrgList` file for each cell:  

``` bash
R

library(devtools)
load_all("bambu")

reads <- "data/demultiplexedONT.bam"
annotations <- "data/Homo_sapiens.GRCh38.91.sorted.gtf"
genome <- "data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
annotations <- prepareAnnotations(annotations)

# Transcript discovery and generate readGrgList for each cell
bambu(reads = reads, annotations = annotations, genome = genome, quant = FALSE, demultiplexed = TRUE,NDR = 1,  
      yieldSize = 1000000)

q()
```

Finally, we can run bambu again for transcript quantification: 

``` bash
### Transcript Quantification for each cell

cellBarcodeDirName <- paste(system.file("extdata", package = "bambu"), "CB", sep = "/")
readGrgListFile <- paste0(cellBarcodeDirName, "/", list.files(cellBarcodeDirName))
rcOutDir <- paste(system.file("extdata", package = "bambu"), "CBReadClass", sep = "/")

annotations <- readRDS("extended_annotations.rds")

se <- bambu(reads = reads, annotations = annotations, genome = genome, ncore = 4, discovery = FALSE, rcOutDir = rcOutDir, 
            readGrgListFile = readGrgListFile)

writeBambuOutput(se, path = "output")

q()
```
The files in the `output` folder is described below:

| Output file name                | Description                                                             |
|:----------------------------|:------------------------------------------|
| extended_annotations.gtf        | Extended transcript & gene annotations for the genome using long reads data.        |
| counts_transcript.txt           | Total read counts estimates for each transcript in each sample (sparse matrix format).        |
| CPM_transcript.txt              | Counts per million (CPM) estimates for each transcript in each sample (sparse matrix format). |
| fullLengthCounts_transcript.txt | Full length read counts estimates for each transcript in each sample (sparse matrix format).  |
| uniqueCounts_transcript.txt                | Unique read counts estimates for each transcript in each sample (sparse matrix format).       |
| txNameToGeneIDMap.txt                 | Gene ID associated to each transcript arranged as in the transcript count estimates          |
| counts_gene.txt                 | Gene read counts estimates for each transcript in each sample.         |
| GeneIDMap.txt                 | Gene ID arranged as in the gene count estimates          |
