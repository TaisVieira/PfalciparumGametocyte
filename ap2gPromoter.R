########################################
## Checks if promoter region of genes ##
##       of interest have AP2-G       ##
##        and/or AP2-G2 motifs        ## 
########################################

library(stringr)

# get promoter region sequence for genes of interest
troph_prog_genes <- read.table("../troph_prog_genes.txt")
gff.path <- "../PlasmoDB-59_Pfalciparum3D7.gff"
gff <- read.table(gff.path, header = F, sep = '\t', quote = NULL)
gff <- gff[which(gff$V3 == 'protein_coding_gene' | gff$V3 == 'pseudogene'), ]
tmp <- strsplit(gff$V9, split = ';')
tmp <- lapply(tmp, function(x){sapply(strsplit(x, split = "="), '[[', 2)})
gene_id <- lapply(tmp, function(x){x[1]})
gene_id <- unlist(gene_id)
ind <- unlist(lapply(troph_prog_genes$V1, function(x){grep(x, gff$V9)}))
genes.int.gff <- gff[ind, ]

genes.int.gff$V10 <- ifelse(genes.int.gff$V7 == '+', genes.int.gff$V4-1000, genes.int.gff$V5+1000) 
genes.int.gff$V11 <- ifelse(genes.int.gff$V7 == '+', genes.int.gff$V4, genes.int.gff$V5)

genes.int.gff$V12 <- ifelse(genes.int.gff$V7 == '+', genes.int.gff$V10, genes.int.gff$V11)
genes.int.gff$V13 <- ifelse(genes.int.gff$V7 == '+', genes.int.gff$V11, genes.int.gff$V10)

chr <- as.factor(genes.int.gff$V1)
tmp <- strsplit(genes.int.gff$V9, split = ';')
tmp <- lapply(tmp, function(x){sapply(strsplit(x, split = "="), '[[', 2)})
gene_id <- unlist(lapply(tmp, function(x){x[1]}))
gene_name <- unlist(lapply(tmp, function(x){x[2]}))
gene_desc <- unlist(lapply(tmp, function(x){x[3]}))
prm_start <- genes.int.gff$V12
prm_end <- genes.int.gff$V13
bed_file <- data.frame(chr, prm_start, prm_end)
bed_file <- bed_file[order(chr), ]
df <- data.frame(chr, gene_id, gene_name, gene_desc, prm_start, prm_end)
df <- df[order(chr), ]

write.csv(df, "../promoters.csv")
write.table(bed_file, "../promoters.bed", row.names = F, 
            col.names = F, sep = '\t', quote = F)

# in comand line: bedtools getfasta -fi genome_file.fasta -bed file.bed -fo file.fasta
# run resulting fasta file in motif finder (BaMM with z-score threshold of 5)
# create an info.file with genes of interest's Gene ID, Genomic Location and Product Description with PlasmoDB

# Inputs: info.file.csv, motif.occurence and promoters.csv
# Output: updated info.file with motifs occurrence
# run once for each motif

info.path<-"../genes.info.csv"
motif.path<-"../promoters_Motif_2/promoters_motif_2.occurrence" #change for each motif
gene.path<-"../promoters.csv"

info.file<-read.csv(info.path)
motif.file<-read.delim(motif.path)
gene.file<-read.csv(gene.path)

gene.file$location<-paste(gene.file$prm_start, gene.file$prm_end, sep='-')
motif.file$location<-sapply(strsplit(motif.file$seq, split=":"),'[[',2)
motif.df<-data.frame()
for (loc in gene.file$location){
  ind<-grep(loc, motif.file$location)
  if (length(ind) != 0){
    motifs<-motif.file$pattern[ind]
    motifs<-paste(motifs, collapse=';')
    ind2<-grep(loc, gene.file$location)
    gene<-gene.file$gene_id[ind2]
    if(ncol(motif.df) == 0){
      motif.df<-data.frame(gene, motifs)
    } else{
      new.row<-data.frame(gene, motifs)
      motif.df = rbind(motif.df, new.row)
    }
  }
}

ind<-match(motif.df$gene, info.file$Gene.ID)
info.file[ind,'AP2.G.motif']<-motif.df$motifs #change column name for each motif

write.csv(info.file, "../genes.info.csv")
