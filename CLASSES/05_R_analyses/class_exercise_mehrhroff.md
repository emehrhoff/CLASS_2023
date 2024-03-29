Class Exercise
================
Erika Mehrhoff
3/17/2023

## Load the libraries you need

## Load functions you need “my\_class\_functions”

## load in your peak files for each replicate of each protein

## Here I am starting to analyze my data for my proteins of interest:

## protein CEBPZ, CHD2, CTCF, ELF1, EP300

## First I will read in each replicate file

``` r
# setting file paths
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/erme3555"
peak_path <- "group/results/bwa/mergedLibrary/macs/broadPeak"
broadpeakfilepath <- file.path(basepath, peak_path)

# import peaks
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# printing out a table of the number of peaks in each file:
# list of how many peaks are in each file before we create consensus peaks.
peak_num <- sapply(peak_list, length) %>% as.data.frame()
# label column
names(peak_num) <- c("num_peaks")
# make dbp name a col.
peak_num <- peak_num %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")

#print summary of peaks
kableExtra::kbl(peak_num, align="c") %>% kableExtra::kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:center;">
dbp
</th>
<th style="text-align:center;">
replicate
</th>
<th style="text-align:center;">
num\_peaks
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center;">
CEBPZ
</td>
<td style="text-align:center;">
R1
</td>
<td style="text-align:center;">
273
</td>
</tr>
<tr>
<td style="text-align:center;">
CEBPZ
</td>
<td style="text-align:center;">
R2
</td>
<td style="text-align:center;">
402
</td>
</tr>
<tr>
<td style="text-align:center;">
CHD2
</td>
<td style="text-align:center;">
R1
</td>
<td style="text-align:center;">
22119
</td>
</tr>
<tr>
<td style="text-align:center;">
CHD2
</td>
<td style="text-align:center;">
R2
</td>
<td style="text-align:center;">
13012
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R1
</td>
<td style="text-align:center;">
36870
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R2
</td>
<td style="text-align:center;">
65526
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R3
</td>
<td style="text-align:center;">
83141
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R4
</td>
<td style="text-align:center;">
55484
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R5
</td>
<td style="text-align:center;">
52373
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R6
</td>
<td style="text-align:center;">
44241
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R7
</td>
<td style="text-align:center;">
71705
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
R8
</td>
<td style="text-align:center;">
68383
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
R1
</td>
<td style="text-align:center;">
25895
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
R2
</td>
<td style="text-align:center;">
35371
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
R3
</td>
<td style="text-align:center;">
23188
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
R4
</td>
<td style="text-align:center;">
15436
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
R1
</td>
<td style="text-align:center;">
38462
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
R2
</td>
<td style="text-align:center;">
59348
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
R3
</td>
<td style="text-align:center;">
4471
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
R4
</td>
<td style="text-align:center;">
6888
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
R5
</td>
<td style="text-align:center;">
10270
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
R6
</td>
<td style="text-align:center;">
882
</td>
</tr>
</tbody>
</table>

## Now I am going to create consensus peaks for each protein

``` r
# list of unique dbps as input parameter
dbps <- unique(sapply(names(peak_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))

# run our function consensus_from_reduced
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)
names(consensus_list) <- dbps

# export consensus peaks to results folder
# setting file path to export
consensus_path <- "CLASS_2023/CLASSES/05_R_analyses/results/"
exportpath <- file.path(basepath, consensus_path)

# export each as .bed file
for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0(exportpath, names(consensus_list)[i], "_consensus_peaks.bed") )}

# clean up consensus and export
# file list:
consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/erme3555/CLASS_2023/CLASSES/05_R_analyses/results/", full.names = T, pattern = "*consensus_peaks.bed")

# lapply (for loop) across consensus file list to add colnames
# The actual col names for .broadPeak are: chr, start, end, name, score, strand
peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))
names(peaks) <- dbps
# double check order by looking at consensus_file_list is same order as dbps

# make chromosomes of interest object
canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")

# let's use lapply with filter funciton to cannonical_chr
peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))

# export to consensus peaks folder
new_filenames <- paste0("results/consensus_peaks/", names(peaks), "_consensus.bed")

for(i in 1:length(peaks)) {
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE) # t for tab separated 
}
```

## Now I am going to make my consensus peaks compatable with UCSC genome browser

``` r
# print out consensus peak files in a results/UCSC directory
# export using naming for ucsc
headers <- paste0("track type=bedGraph name=", names(peaks)) 
new_filenames <- paste0("results/ucsc_directory/", names(peaks), ".bed")

for(i in 1:length(peaks)) {
  # Write the header line
  writeLines(headers[[i]], new_filenames[[i]])
  # Append the broadPeak table data
  
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}
```

## I am curious if my proteins are transcription factors so I will use the annotations in a cell paper I found and see

``` r
# import consensus bed files
consensus_path <- "CLASS_2023/CLASSES/05_R_analyses/results/consensus_peaks"
consensusPeakPath <- file.path(basepath, consensus_path)

consensus_peaks_files <- list.files(consensusPeakPath, 
                                             pattern = "*consensus.bed",
                                             full.names = TRUE)

# lapply with import function to make a list of GRanges
consensus_peaks <- lapply(consensus_peaks_files, rtracklayer::import)

names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/CLASS_2023/erme3555/CLASS_2023/CLASSES/05_R_analyses/results/consensus_peaks/|_consensus.bed","", consensus_peaks_files)

# import gencode Granges
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")


# loading in the number of peaks each DBP has -- using length.
num_peaks_df <- data.frame("dbp" = names(consensus_peaks),
                           "num_peaks" = sapply(consensus_peaks, length))

# total amount of the genome covered by all the peaks for a given DBP.
num_peaks_df$total_peak_length <- sapply(consensus_peaks, function(x) sum(width(x)))

# downloading TF annotation data set
url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"

destination_for_url <- "results/TF_annotations.xlsx"

# to download we can use download.file
download.file(url, destination_for_url)

# import excel sheet to use
human_tfs <- readxl::read_excel("results/TF_annotations.xlsx",
                                sheet = 2, skip = 1)

# rename the 4th column to indicate if it is a TF.
names(human_tfs)[4] <- "is_tf"

# now let's intersect gene names that are in our ChIP data and has TF identity.
length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name))) # 5
```

    ## [1] 5

``` r
# first let's filter and grab the first 4 columns that match DBPs in num_peaks_df
human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]
# which of 4 mapped to downloaded file 

# adding new column names
names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")

# merge
num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T) 

# Let's check how many NAs -- we should have some missing values. 
dim(num_peaks_df[is.na(num_peaks_df$tf),])
```

    ## [1] 0 6

``` r
# if you leave the object name you just created in the environment
# it will print out in the knit
kableExtra::kbl(num_peaks_df, align="c") %>% kableExtra::kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:center;">
dbp
</th>
<th style="text-align:center;">
num\_peaks
</th>
<th style="text-align:center;">
total\_peak\_length
</th>
<th style="text-align:center;">
ensembl\_id
</th>
<th style="text-align:center;">
dbd
</th>
<th style="text-align:center;">
tf
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center;">
CEBPZ
</td>
<td style="text-align:center;">
172
</td>
<td style="text-align:center;">
67337
</td>
<td style="text-align:center;">
ENSG00000115816
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
Yes
</td>
</tr>
<tr>
<td style="text-align:center;">
CHD2
</td>
<td style="text-align:center;">
8815
</td>
<td style="text-align:center;">
4519174
</td>
<td style="text-align:center;">
ENSG00000173575
</td>
<td style="text-align:center;">
Myb/SANT
</td>
<td style="text-align:center;">
No
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
28799
</td>
<td style="text-align:center;">
18606505
</td>
<td style="text-align:center;">
ENSG00000102974
</td>
<td style="text-align:center;">
C2H2 ZF
</td>
<td style="text-align:center;">
Yes
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
11874
</td>
<td style="text-align:center;">
7572457
</td>
<td style="text-align:center;">
ENSG00000120690
</td>
<td style="text-align:center;">
Ets
</td>
<td style="text-align:center;">
Yes
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
568
</td>
<td style="text-align:center;">
500333
</td>
<td style="text-align:center;">
ENSG00000100393
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
No
</td>
</tr>
</tbody>
</table>

## Now I want to compare a protein with a previous analysis

``` r
# go to UCSC genome browser and load in a peak file for a given protein
# load in the data for the same protein from the previous analysis
# compare how your consensus peaks are similar or different to previous analyses

# It appears most of the consensus peaks from this analysis match those found in previous analyses. Scrolling across the genome on UCSC genome browser and looking at the peaks it is evident that almost all of them have similar called peaks.
```

## Now I am going to determine how my peaks for each protein overlap annotations of the genome

## First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

``` r
# pulling out annotations
# all genes
gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 

# exporting all genes file (we will save all the .Rdata too at the end)
rtracklayer::export(gencode_genes, "results/gene_annotations/gencode_genes.gtf")

# mRNA genes (called "protein_coding") in this version of gencode changes sometimes !
mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 

rtracklayer::export(mrna_genes, "results/gene_annotations/mrna_genes.gtf")

# now doing a second index for lncRNA:
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 

rtracklayer::export(lncrna_genes, "results/gene_annotations/lncrna_genes.gtf")

# both mRNA and lncRNA annotations together.
mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
rtracklayer::export(mrna_lncrna_genes, "results/gene_annotations/mrna_lncrna_genes.gtf")

# starting annotation file that we will use moving forward.
lncrna_mrna_genes <- rtracklayer::import("results/gene_annotations/mrna_lncrna_genes.gtf")

# as a dataframe
lncrna_mrna_genes_df <- lncrna_mrna_genes %>% as.data.frame()

lncrna_mrna_promoters <- promoters(lncrna_mrna_genes, upstream = 1000, downstream = 1000)

rtracklayer::export(lncrna_mrna_promoters, "results/gene_annotations/lncrna_mrna_promoters.gtf")

# counting promoter overlaps
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, type = "counts")

# row sum for each DBP
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

# percent promoter overlap
num_peaks_df$percent_overlap <- round((num_peaks_df$peaks_overlapping_promoters/num_peaks_df$num_peaks*100),0)

# find overlaps of promoters for each protein
kableExtra::kbl(num_peaks_df, align="c") %>% kableExtra::kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:center;">
dbp
</th>
<th style="text-align:center;">
num\_peaks
</th>
<th style="text-align:center;">
total\_peak\_length
</th>
<th style="text-align:center;">
ensembl\_id
</th>
<th style="text-align:center;">
dbd
</th>
<th style="text-align:center;">
tf
</th>
<th style="text-align:center;">
peaks\_overlapping\_promoters
</th>
<th style="text-align:center;">
percent\_overlap
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center;">
CEBPZ
</td>
<td style="text-align:center;">
172
</td>
<td style="text-align:center;">
67337
</td>
<td style="text-align:center;">
ENSG00000115816
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
213
</td>
<td style="text-align:center;">
124
</td>
</tr>
<tr>
<td style="text-align:center;">
CHD2
</td>
<td style="text-align:center;">
8815
</td>
<td style="text-align:center;">
4519174
</td>
<td style="text-align:center;">
ENSG00000173575
</td>
<td style="text-align:center;">
Myb/SANT
</td>
<td style="text-align:center;">
No
</td>
<td style="text-align:center;">
8088
</td>
<td style="text-align:center;">
92
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
28799
</td>
<td style="text-align:center;">
18606505
</td>
<td style="text-align:center;">
ENSG00000102974
</td>
<td style="text-align:center;">
C2H2 ZF
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
3848
</td>
<td style="text-align:center;">
13
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
11874
</td>
<td style="text-align:center;">
7572457
</td>
<td style="text-align:center;">
ENSG00000120690
</td>
<td style="text-align:center;">
Ets
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
9304
</td>
<td style="text-align:center;">
78
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
568
</td>
<td style="text-align:center;">
500333
</td>
<td style="text-align:center;">
ENSG00000100393
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
No
</td>
<td style="text-align:center;">
68
</td>
<td style="text-align:center;">
12
</td>
</tr>
</tbody>
</table>

## Results:

## 1) What can you determine from these overlaps?

### CEBPZ binds at promoters a lot (124%)

### CHD2 binds at promoters a lot (92%)

### CTCF does not bind at promoters a lot (13%)

### ELF1 binds at promoters a lot (78%)

### EP300 does not bind at promoters a lot (12%)

## Now I want to compare the overlaps with lncRNA and mRNA promoters seperately

``` r
# last handy annotation will be lncRNA and mRNA gene IDs to subset
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]

# break these promoters into two groups "lncrna" and "mrna" using the gene_id objects we made above to index and separate them.
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

kableExtra::kbl(num_peaks_df[,c(1,2,5,6,7,9,10)], align="c") %>% kableExtra::kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:center;">
dbp
</th>
<th style="text-align:center;">
num\_peaks
</th>
<th style="text-align:center;">
dbd
</th>
<th style="text-align:center;">
tf
</th>
<th style="text-align:center;">
peaks\_overlapping\_promoters
</th>
<th style="text-align:center;">
peaks\_overlapping\_lncrna\_promoters
</th>
<th style="text-align:center;">
peaks\_overlapping\_mrna\_promoters
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center;">
CEBPZ
</td>
<td style="text-align:center;">
172
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
213
</td>
<td style="text-align:center;">
34
</td>
<td style="text-align:center;">
179
</td>
</tr>
<tr>
<td style="text-align:center;">
CHD2
</td>
<td style="text-align:center;">
8815
</td>
<td style="text-align:center;">
Myb/SANT
</td>
<td style="text-align:center;">
No
</td>
<td style="text-align:center;">
8088
</td>
<td style="text-align:center;">
1694
</td>
<td style="text-align:center;">
6394
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
28799
</td>
<td style="text-align:center;">
C2H2 ZF
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
3848
</td>
<td style="text-align:center;">
1337
</td>
<td style="text-align:center;">
2511
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
11874
</td>
<td style="text-align:center;">
Ets
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
9304
</td>
<td style="text-align:center;">
1973
</td>
<td style="text-align:center;">
7331
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
568
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
No
</td>
<td style="text-align:center;">
68
</td>
<td style="text-align:center;">
39
</td>
<td style="text-align:center;">
29
</td>
</tr>
</tbody>
</table>

## Results:

## 1) What is the difference in overlaps between mRNA and lncRNA promoters

### All but 1 of the 5 proteins (EP300) bind promoters more in mRNA than lncRNA.

## Now I am going to test if there is more binding over gene bodies than promoters

## I will seperate lncRNA and mRNA gene bodies to find the overlaps

``` r
# Finding overlaps with gene_bodies
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                consensus_peaks, 
                                                type = "counts")

# extract the overlaps the same way as promoters above
# All gene bodies
num_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)

# lncRNA gene bodies 
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
num_peaks_df$peaks_overlapping_mrna_genebody <- 
  rowSums(genebody_peak_counts[,mrna_gene_ids])

# export results
write_csv(num_peaks_df, "results/num_peaks_df.csv")

kableExtra::kbl(num_peaks_df[,c(1:7,11)], align="c") %>% kableExtra::kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:center;">
dbp
</th>
<th style="text-align:center;">
num\_peaks
</th>
<th style="text-align:center;">
total\_peak\_length
</th>
<th style="text-align:center;">
ensembl\_id
</th>
<th style="text-align:center;">
dbd
</th>
<th style="text-align:center;">
tf
</th>
<th style="text-align:center;">
peaks\_overlapping\_promoters
</th>
<th style="text-align:center;">
peaks\_overlapping\_genebody
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center;">
CEBPZ
</td>
<td style="text-align:center;">
172
</td>
<td style="text-align:center;">
67337
</td>
<td style="text-align:center;">
ENSG00000115816
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
213
</td>
<td style="text-align:center;">
207
</td>
</tr>
<tr>
<td style="text-align:center;">
CHD2
</td>
<td style="text-align:center;">
8815
</td>
<td style="text-align:center;">
4519174
</td>
<td style="text-align:center;">
ENSG00000173575
</td>
<td style="text-align:center;">
Myb/SANT
</td>
<td style="text-align:center;">
No
</td>
<td style="text-align:center;">
8088
</td>
<td style="text-align:center;">
10453
</td>
</tr>
<tr>
<td style="text-align:center;">
CTCF
</td>
<td style="text-align:center;">
28799
</td>
<td style="text-align:center;">
18606505
</td>
<td style="text-align:center;">
ENSG00000102974
</td>
<td style="text-align:center;">
C2H2 ZF
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
3848
</td>
<td style="text-align:center;">
22208
</td>
</tr>
<tr>
<td style="text-align:center;">
ELF1
</td>
<td style="text-align:center;">
11874
</td>
<td style="text-align:center;">
7572457
</td>
<td style="text-align:center;">
ENSG00000120690
</td>
<td style="text-align:center;">
Ets
</td>
<td style="text-align:center;">
Yes
</td>
<td style="text-align:center;">
9304
</td>
<td style="text-align:center;">
13817
</td>
</tr>
<tr>
<td style="text-align:center;">
EP300
</td>
<td style="text-align:center;">
568
</td>
<td style="text-align:center;">
500333
</td>
<td style="text-align:center;">
ENSG00000100393
</td>
<td style="text-align:center;">
Unknown
</td>
<td style="text-align:center;">
No
</td>
<td style="text-align:center;">
68
</td>
<td style="text-align:center;">
439
</td>
</tr>
</tbody>
</table>

## Results:

## 1) Do my proteins have more overlaps with promoters or genebodies?

### All but 1 protein (CEBPZ) overlaps more with genebodies than promoters.

## It is nice and all to find overlaps, but I am interested in how many proteins bind a specific promoter. I will use my handy “occurence” parameter in " count peaks per feature"

``` r
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, 
                                               type = "occurrence") 

# check that all lncrna & mrna genes are accounted for:
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

# will use in future
write.table(promoter_peak_occurence, "results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# make into dataframe
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))

# now call in one csv
write_csv(peak_occurence_df, "results/peak_occurence_dataframe.csv")

# summary of number of dbp binding promoter
sum_tab <- table(peak_occurence_df$number_of_dbp)
names(dimnames(sum_tab)) <- c("Max number of proteins bound specific promoter")
sum_tab
```

    ## Max number of proteins bound specific promoter
    ##     0     1     2     3     4 
    ## 24271  5895  5363  1270    15

## Results: I find the max number of proteins on a promoter to be 4.

## Now I want to start plotting my results

## First I will see if there is a realtionship between peak number and total DNA covered

``` r
ggplot(num_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length)) +
  geom_point() +
  labs(x="Number of Peaks",
       y="Total Peak Length",
       title= "Number of Peaks and the Total DNA Covered") +
  geom_smooth(method = "lm", se=F)
```

![](class_exercise_mehrhroff_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Now I want to color my plot by whether the protein is a TF or not.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length, color=tf)) +
  geom_point() +
  labs(x="Number of Peaks",
       y="Total Peak Length",
       title= "Number of Peaks and the Total DNA Covered") 
```

![](class_exercise_mehrhroff_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## I want to make a histogram of the number of peaks for each of my proteins

``` r
ggplot(num_peaks_df, aes(x = num_peaks)) +
  geom_histogram(bins=5, fill="#6495ED", color="black") +
  labs(x="Number of Peaks",
       y="Frequency",
       title= "Number of Peaks for all the Proteins")
```

![](class_exercise_mehrhroff_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Now I want to facet this by the type of DNA binding domain my protein has.

``` r
#ggplot(num_peaks_df, aes(x = num_peaks)) +
  #geom_histogram(bins=7) +
  #facet_wrap(~dbd)
# If you take previous graph and facet wrap it by dbp it won't look right because each histogram would only have 1 piece of information.
# Could represent as histogram that has colored bars for each dbd or a bar plot

ggplot(num_peaks_df, aes(x = num_peaks, fill=dbd)) +
  geom_histogram(bins=5) +
  labs(x="Number of Peaks",
       y="Frequency",
       title= "Number of Peaks for all the Proteins")
```

![](class_exercise_mehrhroff_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggplot(num_peaks_df, aes(
      x = dbd,
      y = num_peaks, fill=dbd)) +
        geom_bar(stat = "identity") +
  labs(x="Type of DNA Binding Protein",
       y="Number of Peaks",
       title= "Number of Peaks by Type of DNA Binding Protein") 
```

![](class_exercise_mehrhroff_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

## Cool now I am ready to send my result to my collaborator as a Knitted document
