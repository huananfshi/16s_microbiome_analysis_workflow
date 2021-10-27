library(dada2)
library(biomformat)
library(ALDEx2)
library(DESeq2)
library(phyloseq)
library(phangorn)
library(vegan)
library(DECIPHER)
dada2_workflow_16s=function(path,site){
  print(path)
  print(site)
  list.files(paste(path,site,sep=""))
  fnFs=sort(list.files(paste(path,site,sep=""), pattern = ".1.fq.bz2",full.names = TRUE))
  fnRs=sort(list.files(paste(path,site,sep=""), pattern = ".2.fq.bz2",full.names = TRUE))
  #16 EBC sample 1 only contain one end of the pair end...
  sample.names = sapply(strsplit(basename(fnFs),".1.fq"),`[`,1)
  plotQualityProfile(fnFs[1:2])
  plotQualityProfile(fnRs[1:2])
  filtFs = file.path(paste(path,site,sep=""),"filtered",paste0(sample.names,"_F_filt.fastq.gz"))
  filtRs = file.path(paste(path,site,sep=""),"filtered",paste0(sample.names,"_R_filt.fastq.gz"))
  names(filtFs)=sample.names
  names(filtRs)=sample.names
  out = filterAndTrim(fnFs,filtFs,fnRs,filtRs,truncLen = c(240,160),maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,compress = T,multithread = T)
  errF = learnErrors(filtFs,multithread = T)
  errR = learnErrors(filtRs,multithread = T)
  dadaFs = dada(filtFs,err = errF, multithread = T)
  dadaRs = dada(filtRs,err = errR, multithread = T)
  mergers=mergePairs(dadaFs,filtFs,dadaRs,filtRs,verbose = T)
  seqtab = makeSequenceTable(mergers)
  head(seqtab)
  seqtab.nochim = removeBimeraDenovo(seqtab, method = "consensus", multithread = T, verbose = T)
  sum(seqtab.nochim)/sum(seqtab)
  asv_seqs = colnames(seqtab.nochim)
  asv_headers = vector(dim(seqtab.nochim)[2],mode = "character")
  for (i in 1: dim(seqtab.nochim)[2]){
    asv_headers[i]=paste(">ASV",i,sep = "_")
  }
  asv_fasta = c(rbind(asv_headers,asv_seqs))
  str(asv_fasta)
  write(asv_fasta,paste(path,"dada2_output/asv_",site,".fna",sep=""))
  asv_tab = t(seqtab.nochim)
  rownames(asv_tab) = sub(">","",asv_headers)
  st.biom = make_biom(asv_tab)
  write_biom(st.biom,paste(path,"dada2_output/asv_",site,".biom",sep=""))
  
  #taxonomy assignment
  taxa =assignTaxonomy(seqtab.nochim,"~/Desktop/16s_functional_profiles/dada2_output/tax_ref/silva_nr_v138_train_set.fa.gz",multithread = T)
  
  taxa_frame = data.frame(taxa)
  rownames(taxa_frame) = asv_headers
  # taxa_species = assignSpecies(seqtab.nochim, "~/Desktop/16s_functional_profiles/dada2_output/tax_ref/silva_species_assignment_v138.fa.gz")
  # rownames(taxa_species) = asv_headers
  # taxa_species_frame = data.frame(taxa_species)
  # taxa_frame$Species = taxa_species_frame$Species
  
  # write.table(taxa_frame,"~/Box Sync/EODF_analysis/intermediate_files/16s/dada2_output/otu_taxonomy_all.tsv",sep = "\t",quote=F,row.names = T)
  # taxa_frame$taxonomy = paste(taxa_frame$Kingdom,taxa_frame$Phylum,taxa_frame$Class,taxa_frame$Order,taxa_frame$Family,taxa_frame$Genus,taxa_frame$Species,sep=";")
  # taxa_frame$taxonomy = gsub(";NA","",taxa_frame$taxonomy)
  # asv_tab_taxa = cbind(asv_tab,taxa_frame$taxonomy)
  # head(asv_tab_taxa)
  # write_biom(st.biom,"~/Box Sync/EODF_analysis/intermediate_files/16s/dada2_output/otu_all_taxa.biom")
  
  taxa_frame$taxonomy = paste(taxa_frame$Kingdom,taxa_frame$Phylum,taxa_frame$Class,taxa_frame$Order,taxa_frame$Family,taxa_frame$Genus,sep=";")
  taxa_frame$taxonomy = gsub("(;NA)+",";unclassified",taxa_frame$taxonomy)
  rownames(taxa_frame)=gsub(">","",rownames(taxa_frame))
  write.csv(taxa_frame,paste(path,"dada2_output/asv_genus_",site,".csv",sep=""),quote=F,row.names = T) #write asv_genus list
  return_list = list("seqtabNoC" = seqtab.nochim, "taxTab"= taxa, "asv_tab"= asv_tab)
  return(return_list)
  }
