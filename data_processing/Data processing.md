## NADs identification

#### 1. Sort and index aligned bam files

##### After aligned the sequencing data of WGS, NS, nHi-C and *in situ* Hi-C, we need to sort and index the aligned bam files and be ready for calculating the coverage. 

```shell
while read -a line do samtools sort $line > ${line%%.*}.sort.bam & done < $1

while read -a line do samtools index $line & > done < $1
```

#### 2. Transform bam files to [bedgraph](http://genome.ucsc.edu/goldenpath/help/bedgraph.html) files

##### Then spliting the hg19 genome into 1kb bins and we can use samtools bedcov function to calculate coverage of each genome bins.

```shell
while read -a line do samtools bedcov -Q 30 hg19_1k.bed $line > ${line%%.*}.q30.bedgraph & done < bam_list
```

#### 3. Bedgraph files normalization

##### To correct the coverage bias caused by sequencing depth, we normalized the bin coverage by divide the total number of aligned reads.

```R
library(data.table)

files <- list.files("./bedgraph/", pattern = "*q30.1k", full.names = T)

for (i  in 1:length(files)) {

	mat <- fread(files[i])

	mat$V6 <- as.numeric(mat$V6)/sum(as.numeric(mat$V6)) * 100000000

	sample <- basename(files[i])

	write.table(mat, paste0("./bedgraph/", sample, ".norm.BedGraph"),

					row.names = F, col.names = F, quote = F, sep = "\t")

}
```

#### 4. Calculate log2nHi-C/Hi-C and log2NS/WGS ratio

##### To identify NADs with HMM model, we need to prepare the file of log2nHi-C/Hi-C and log2NS/WGS ratio, the Rsciprt for calculating log2NS/WGS ratio are showed here:

```R
setwd("./dnabam/bedgraph/")

wgs1 <- fread("hela_wgs1.q30.1k.bedgraph.norm.BedGraph")

wgs2 <- fread("hela_wgs2.q30.1k.bedgraph.norm.BedGraph")

n1 <- fread("hela_n1_wgs.q30.1k.bedgraph.norm.BedGraph")

n2 <- fread("hela_n2_wgs.q30.1k.bedgraph.norm.BedGraph")

wgs <- cbind(wgs1[,c(1:3,4)], wgs2[,4])

wgs$mean <- rowMeans(wgs[,3:4])

mat <- cbind(wgs[,c(1,5)], n2[,4])

names(mat) <- c("chromosome", "position", "wgs","n2")

mat$ratio <- log2(mat$n2/mat$wgs)

loci <- c(which(is.na(mat$ratio)), which(is.infinite(mat$ratio)))

mat$ratio[loci] <- mat$n2[loci] - mat$wgs[loci]

mat <- mat[-which(mat$wgs == 0 & mat$n2 == 0), ]

saveRDS(mat, 'mat.rds')
```

#### 5. HMM calling NADs

##### Now, we can use the HMM model to call NADs. We have defined a function [getNAD()](https://github.com/ChengLiLab/nHi-C/blob/main/data_processing/getNAD.R) to process each log2nHi-C/Hi-C and log2NS/WGS ratio matrix, the Rsciprt for calling NADs using  log2NS/WGS ratio matrix are showed here:

```R
source('getNAD.R')

mat <- readRDS('mat.rds')

wgmat <- split(mat, f = mat$chromosome)

wgregions <- data.frame(chr = as.character(), start = as.numeric(), end = as.numeric())

for(i in 1:23){

	if(i == 23){chr = "chrX"}

	else{chr = paste0("chr", i)}

	chrmat <- wgmat[[chr]]

	wgregions <- rbind(wgregions, getNAD(chrmat))

}

write.table(wgregions, "nad.list", row.names = F, col.names = F, quote = F, sep = "\t")
```

#### 6. Merge and filter NADs

##### Finally, we merged the NADs within 100kb and removed the NADs whose length were   less than 10kb.

```shell
while read ia line do bedtools merge -d 100000 -i $line > merge_$line done < nad.list

while read -a line do awk '$3-$2>=10000' $line.merged.bed > $line.merged.filter.bed done <  nad_merge.list
```
