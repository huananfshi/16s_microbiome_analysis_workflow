#loading required packages
library(dada2)
library(biomformat)
path = "~/.../qiime2_data" # the output folder of qimme2
list.files(path)
fnFs=sort(list.files(path, pattern = "_R1_001.fastq.gz",full.names = TRUE))
fnRs=sort(list.files(path, pattern = "_R2_001.fastq.gz",full.names = TRUE))
sample.names = sapply(strsplit(basename(fnFs),"_"),`[`,1)
#inspect real quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
#Filter and trim
filtFs = file.path(path,"filtered",paste0(sample.names,"_F_filt.fastq.gz"))
filtRs = file.path(path,"filtered",paste0(sample.names,"_R_filt.fastq.gz"))
names(filtFs)=sample.names
names(filtRs)=sample.names
out = filterAndTrim(fnFs,filtFs,fnRs,filtRs,truncLen = c(240,160),maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,compress = T,multithread = T)
#Learn the Error Rates
errF = learnErrors(filtFs,multithread = T)
errR = learnErrors(filtRs,multithread = T)
#Sample Inference
dadaFs = dada(filtFs,err = errF, multithread = T)
dadaRs = dada(filtRs,err = errR, multithread = T)
#Merge paired reads
mergers=mergePairs(dadaFs,filtFs,dadaRs,filtRs,verbose = T)
#Construct sequence table
seqtab = makeSequenceTable(mergers)
head(seqtab)
#Remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method = "consensus", multithread = T, verbose = T)
sum(seqtab.nochim)/sum(seqtab)
#generate ASVs biom file and ASVs sequences (.fasta) for PICRUSt2
asv_seqs = colnames(seqtab.nochim)
asv_headers = vector(dim(seqtab.nochim)[2],mode = "character")
for (i in 1: dim(seqtab.nochim)[2]){
asv_headers[i]=paste(">ASV",i,sep = "_")
}
asv_fasta = c(rbind(asv_headers,asv_seqs))
write(asv_fasta,"~/.../dada2_output/asv.fna")
asv_tab = t(seqtab.nochim)
rownames(asv_tab) = sub(">","",asv_headers)
st.biom = make_biom(asv_tab)
write_biom(st.biom,"~/.../dada2_output/asv.biom")

#taxonomy assignment
taxa =assignTaxonomy(seqtab.nochim,"~/.../tax_ref/silva_nr_v138_train_set.fa.gz",multithread = T) #path for tax_ref databases
rownames(taxa) = asv_headers
taxa_frame = data.frame(taxa)
taxa_species = assignSpecies(seqtab.nochim, "~/...t/tax_ref/silva_species_assignment_v138.fa.gz") #optional, for species level
rownames(taxa_species) = asv_headers
taxa_species_frame = data.frame(taxa_species)
taxa_frame$Species = taxa_species_frame$Species

###write.table(taxa_frame,"~/.../dada2_output/asv_taxonomy.tsv",sep = "\t",quote=F,row.names = T)
###taxa_frame$taxonomy = paste(taxa_frame$Kingdom,taxa_frame$Phylum,taxa_frame$Class,taxa_frame$Order,taxa_frame$Family,taxa_frame$Genus,taxa_frame$Species,sep=";")
###taxa_frame$taxonomy = gsub(";NA","",taxa_frame$taxonomy)
###asv_tab_taxa = cbind(asv_tab,taxa_frame$taxonomy)
###head(asv_tab_taxa)
###write_biom(st.biom,"~/Desktop/16s_functional_profiles/dada2_output/asv_taxa.biom") # write biom files with taxa

taxa_frame$taxonomy = paste(taxa_frame$Kingdom,taxa_frame$Phylum,taxa_frame$Class,taxa_frame$Order,taxa_frame$Family,asv_taxa$Genus,sep=";")
taxa_frame$taxonomy = gsub(";NA",";unclassified",taxa_frame$taxonomy)
rownames(taxa_frame)=gsub(">","",rownames(taxa_frame))
write.csv(taxa_frame,"~/.../dada2_output/asv_genus.csv",quote=F,row.names = T) #write asv_genus list
