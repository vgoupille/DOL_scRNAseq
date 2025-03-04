library(reticulate)
library(dplyr)
library(cowplot)
library(Seurat) ##Note this code was written for a previous version of Seurat (4.1.1)
library(ggplot2)
library(Matrix)
library(forcats)
library(akmedoids)
library(RColorBrewer)

#--yeast data processing--------------------------------------------------------
#function to create Seurat Object for downstream analysis from counts, features, barcodes, and other metadata
create_seurat_object_from_seq_files <- function(data_dir, sublibrary, ribosome_removal, elbowcutoff, doublet_percentage, experiment){
  
  setwd(list.dirs(data_dir))
  
  data <- readMM('UniqueAndMult-Uniform.mtx')
  genes <- read.table('features.tsv')
  barcodes <- read.table('barcodes.tsv')
  genes <- genes$V1
  genes <- convert_gene_names_yeast(genes, TRUE)
  barcodes <- barcodes$V1
  
  rownames(data) <- genes
  colnames(data) <- barcodes
  #data[which(rowSums(data) > 0),]
  
  #use filtered cells
  filteredData <- round_one_bc_collapse(data, elbowcutoff)
  
  #deal with ribosomal reads
  
  if (ribosome_removal == 'rRNA') {
    ribosomes <- read.csv('C:/Users/lbrettne/Desktop/Solo.out/ribosomes.csv')
    colnames(ribosomes)[1] <- 'name1'
    rRNA <- ribosomes[which(ribosomes$type == 'rRNA'),]
    mRNA <- filteredData[which(!(rownames(filteredData) %in% rRNA$name1)),]
    mRNA <- mRNA[which(!(rownames(mRNA) %in% rRNA$name2)),]
    
    #convert to Seurat object
    SO <- CreateSeuratObject(counts = mRNA, project = experiment)
    
  } else {
    
    #convert to Seurat object
    SO <- CreateSeuratObject(counts = filteredData, project = experiment)
    
  }
  
  #add metadata
  if (doublet_percentage == 0) {
    
    wells <- assign_cell_wells(colnames(SO@assays$RNA))
    
    SO$well <- wells$well
    SO$sublibrary <- sublibrary
    SO$experiment <- experiment
    
  } else {
    
    wells <- assign_cell_wells(colnames(SO@assays$RNA))
    species <- assign_cell_species(SO@assays$RNA,doublet_percentage,"rRNA")
    
    celldata <- cbind(wells$well,species$s288c_reads, species$sc5314_reads, species$species)
    colnames(celldata) <- c('well','s288c_reads','sc5314_reads','species')
    celldata <- as.data.frame(celldata)
    row.names(celldata) <- wells$barcode
    
    SO$well <- celldata$well
    SO$s288c_reads <- celldata$s288c_reads
    SO$sc5314_reads <- celldata$sc5314_reads
    SO$species <- celldata$species
    SO$sublibrary <- sublibrary
    SO$experiment <- experiment 
  }
  
  return(SO)
}

#Yeast sequencing sublibraries from the same experiment
sample2 <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GrowthCurve/sample2_reseq/GeneFull/raw/','sample2rs','none',0.85,0,"3")
sample3 <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GrowthCurve/sample3_reseq/GeneFull/raw/','sample3rs','none',0.921,0,"3")

#Same datasets as above with rRNA removed
sample2.nr <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GrowthCurve/sample2_reseq/GeneFull/raw/','sample2rs','rRNA',0.85,0,"3")
sample3.nr <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GrowthCurve/sample3_reseq/GeneFull/raw/','sample3rs','rRNA',0.921,0,"3")

#functions adds cell metadata based on which well barcodes are detected
assign_metadata_yeast <- function(raw.genefull) {         

  s1 <- merge(subset(raw.genefull, subset = well == 'A1'), y = c(subset(raw.genefull, subset = well == 'A2')))
  s1$cond <- 's1'
  s1$hour <- '17.25'
  s1$celltype <- 's288c'
  s2 <- merge(subset(raw.genefull, subset = well == 'A3'), y = c(subset(raw.genefull, subset = well == 'A4')))
  s2$cond <- 's2'
  s2$hour <- '17.25'
  s2$celltype <- 's288c'
  s3 <- merge(subset(raw.genefull, subset = well == 'A5'), y = c(subset(raw.genefull, subset = well == 'A6')))
  s3$cond <- 's3'
  s3$hour <- '17.25'
  s3$celltype <- 'by4741'
  s4 <- merge(subset(raw.genefull, subset = well == 'A7'), y = c(subset(raw.genefull, subset = well == 'A8')))
  s4$cond <- 's4'
  s4$hour <- '17.25'
  s4$celltype <- 'by4741'
  s5 <- merge(subset(raw.genefull, subset = well == 'A9'), y = c(subset(raw.genefull, subset = well == 'A10')))
  s5$cond <- 's5'
  s5$hour <- '19.25'
  s5$celltype <- 's288c'
  s6 <- merge(subset(raw.genefull, subset = well == 'A11'), y = c(subset(raw.genefull, subset = well == 'A12')))
  s6$cond <- 's6'
  s6$hour <- '19.25'
  s6$celltype <- 's288c'
  s7 <- merge(subset(raw.genefull, subset = well == 'B1'), y = c(subset(raw.genefull, subset = well == 'B2')))
  s7$cond <- 's7'
  s7$hour <- '19.25'
  s7$celltype <- 'by4741'
  s8 <- merge(subset(raw.genefull, subset = well == 'B3'), y = c(subset(raw.genefull, subset = well == 'B4')))
  s8$cond <- 's8'
  s8$hour <- '19.25'
  s8$celltype <- 'by4741'
  s9 <- merge(subset(raw.genefull, subset = well == 'B5'), y = c(subset(raw.genefull, subset = well == 'B6')))
  s9$cond <- 's9'
  s9$hour <- '21.25'
  s9$celltype <- 's288c'
  s10 <- merge(subset(raw.genefull, subset = well == 'B7'), y = c(subset(raw.genefull, subset = well == 'B8')))
  s10$cond <- 's10'
  s10$hour <- '21.25'
  s10$celltype <- 's288c'
  s11 <- merge(subset(raw.genefull, subset = well == 'B9'), y = c(subset(raw.genefull, subset = well == 'B10')))
  s11$cond <- 's11'
  s11$hour <- '21.25'
  s11$celltype <- 'by4741'
  s12 <- merge(subset(raw.genefull, subset = well == 'B11'), y = c(subset(raw.genefull, subset = well == 'B12')))
  s12$cond <- 's12'
  s12$hour <- '21.25'
  s12$celltype <- 'by4741'
  s13 <- merge(subset(raw.genefull, subset = well == 'C1'), y = c(subset(raw.genefull, subset = well == 'C2')))
  s13$cond <- 's13'
  s13$hour <- '23.25'
  s13$celltype <- 's288c'
  s14 <- merge(subset(raw.genefull, subset = well == 'C3'), y = c(subset(raw.genefull, subset = well == 'C4')))
  s14$cond <- 's14'
  s14$hour <- '23.25'
  s14$celltype <- 's288c'
  s15 <- merge(subset(raw.genefull, subset = well == 'C5'), y = c(subset(raw.genefull, subset = well == 'C6')))
  s15$cond <- 's15'
  s15$hour <- '23.25'
  s15$celltype <- 'by4741'
  s16 <- merge(subset(raw.genefull, subset = well == 'C7'), y = c(subset(raw.genefull, subset = well == 'C8')))
  s16$cond <- 's16'
  s16$hour <- '23.25'
  s16$celltype <- 'by4741'
  
  growcurv <- merge(s1,y = c(s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16))

  return(growcurv)
}

growcurvy <- assign_metadata_yeast(merge(sample2, y = sample3))
growcurvy.nr <- assign_metadata_yeast(merge(sample2.nr, y = sample3.nr))

#subset data into individual timepoints based on approximate sample time (9AM, 11AM, 1PM, 3PM)
nine <- subset(growcurvy, subset = hour == '17.25')
eleven <- subset(growcurvy, subset = hour == '19.25')
one <- subset(growcurvy, subset = hour == '21.25')
three <- subset(growcurvy, subset = hour == '23.25')

#subset data minus rRNA into individual timepoints
nine.nr <- subset(growcurvy.nr, subset = hour == '17.25')
eleven.nr <- subset(growcurvy.nr, subset = hour == '19.25')
one.nr <- subset(growcurvy.nr, subset = hour == '21.25')
three.nr <- subset(growcurvy.nr, subset = hour == '23.25')

#normalize and scale mRNA datasets
scale_normalize_seurat_data <- function(seuratobject, nvariableFeatures, npcs, features){
  
  integrated <- NormalizeData(seuratobject, verbose = FALSE)
  integrated <- FindVariableFeatures(integrated, selection.method = "vst", 
                                     nfeatures = nvariableFeatures, verbose = FALSE)
  
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE, features = features)
  #integrated <- JackStraw(integrated)
  #integrated <- ScoreJackStraw(integrated, dims = 1:npcs)
  #p <- plot_grid(ElbowPlot(integrated), JackStrawPlot(integrated, dims = 1:npcs))
  p <- ElbowPlot(integrated)
  print(p)
  return(integrated)
}

nine.nr.int <- scale_normalize_seurat_data(nine.nr, 3000, 20, NULL)
eleven.nr.int <- scale_normalize_seurat_data(eleven.nr, 3000, 20, NULL)
one.nr.int <- scale_normalize_seurat_data(one.nr, 3000, 20, NULL)
three.nr.int <- scale_normalize_seurat_data(three.nr, 3000, 20, NULL)

#calculate counts and statistics for RNA types, including rRNA, mRNA, etc.
assign_RNA_types_y <- function(SO, SO2, group_name, timepoint) {
  
  ribosomes <- read.csv('C:/Users/lbrettne/Desktop/Solo.out/ribosomes.csv')
  names(ribosomes)[1] = "name1"
  rRNA.genes <- ribosomes[which(ribosomes$type == 'rRNA'),]
  colnames(rRNA.genes) <- c("name1", "name2", "type", "location")
  tRNA.genes <- ribosomes[which(ribosomes$type == 'tRNA'),]
  colnames(tRNA.genes) <- c("name1", "name2", "type", "location")
  rProt.genes <- ribosomes[which(ribosomes$type == 'rProtein'),]
  colnames(rProt.genes) <- c("name1", "name2", "type", "location")
  
  rRNA  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rRNA.genes$name1 | rownames(SO@assays$RNA) %in% rRNA.genes$name2),])
  tRNA  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% tRNA.genes$name1 | rownames(SO@assays$RNA) %in% tRNA.genes$name2),])
  mRNA  <- SO@assays$RNA[which(!(rownames(SO@assays$RNA) %in% rRNA.genes$name1 | rownames(SO@assays$RNA) %in% rRNA.genes$name2)),]
  mRNA  <- mRNA <- mRNA[which(!(rownames(mRNA) %in% tRNA.genes$name1 | rownames(mRNA) %in% tRNA.genes$name2)),]
  mRNA  <- colSums(mRNA)
  rProt <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rProt.genes$name1 | rownames(SO@assays$RNA) %in% rProt.genes$name2),])
  rProt.ratio <- rProt/mRNA
  totalRNA <- colSums(SO@assays$RNA)
  rRNA.ratio <- rRNA/totalRNA
  bulkrr <- log10(sum(rRNA)/sum(totalRNA)*10^6)
  mRNA.ratio <- mRNA/totalRNA
  mRNA.rRNA <- mRNA/rRNA
  ploidy <- SO$celltype
  
  class <- read.csv('C:/Users/lbrettne/Desktop/RiBi_iESR_CC_genes.csv')
  class$Gene <- convert_gene_names_yeast(class$Gene, FALSE)
  RP.genes <- class[which(class$Group == 'RP'),]
  iESR.genes <- class[which(class$Group == 'iESR'),]
  RiBi.genes <- class[which(class$Group == 'RiBi'),]
  
  RP  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RP.genes$Gene),])
  rpcv <- sd(RP[which(RP>0)])/mean(RP[which(RP>0)])
  iESR  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% iESR.genes$Gene),])
  RiBi  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RiBi.genes$Gene),])
  totalmRNA <- colSums(SO2@assays$RNA)

  bulkrp <- log10(sum(rProt)/sum(totalRNA)*10^6)
  rcv.log <- sd(log(rRNA))/mean(log(rRNA))
  rpcv.log <- sd(log(rProt[which(rProt.ratio >0)]))/mean(log(rProt[which(rProt.ratio >0)]))
  
  fanorr <- sd(rRNA)^2/mean(rRNA)
  fanorp <- sd(rProt[which(rProt>0)])^2/mean(rProt[which(rProt>0)])
  
  rcv <- sd(rRNA)/mean(rRNA)
  mcv <- sd(mRNA)/mean(mRNA)
  rr.med.dev.med <- median(abs(rRNA - median(rRNA)))
  #rpcv <- sd(rProt[which(rProt > 0)])/mean(rProt[which(rProt > 0)])
  rp.med.dev.med <- median(abs(rProt[which(rProt > 0)] - median(rProt[which(rProt > 0)])))

  Gxt <- c( 1.049103, 0.7824108, 0.3477437, 0.09317876)
  
  temp <- data.frame(group_name = group_name, 
                     rRNA = rRNA, 
                     tRNA = tRNA, 
                     mRNA = mRNA, 
                     rProt = rProt,
                     totalRNA = totalRNA,
                     rProt.ratio = rProt.ratio,
                     rRNA.ratio = rRNA.ratio,
                     mRNA.ratio = mRNA.ratio, 
                     mRNA.rRNA = mRNA.rRNA,
                     bulkrr = bulkrr,
                     GpH = Gxt[timepoint],
                     bulkrp = bulkrp,
                     rcv = rcv,
                     mcv = mcv,
                     rpsd = sd(rProt),
                     rpcv = rpcv,
                     rcv.log = rcv.log,
                     rpcv.log = rpcv.log,
                     rm = mean(rRNA),
                     rrm = mean(rRNA.ratio),
                     rsd = sd(rRNA),
                     rrsd = sd(rRNA.ratio),
                     rpm = mean(rProt),
                     rprm = mean(rProt.ratio),
                     rpsd = sd(rProt[which(rProt > 0)]),
                     rprsd = sd(rProt.ratio),
                     rRNA.TPM = (sum(rRNA)/sum(totalRNA))*1e6,
                     rProt.TPM = (sum(rProt)/sum(totalRNA))*1e6,
                     rr.med.dev.med = rr.med.dev.med,
                     rp.med.dev.med = rp.med.dev.med,
                     ploidy = ploidy,
                     gm.rRNA = exp(mean(log(rRNA))),
                     gm.rProt = exp(mean(log(rProt+1))),
                     gm.rr = exp(mean(log(rRNA.ratio))),
                     med.rRNA = median(rRNA),
                     med.rr = median (rRNA.ratio),
                     med.rProt = median(rProt),
                     fanorr = fanorr,
                     fanorp = fanorp) 
  
  return(temp)
}

nine.rna <- assign_RNA_types_y(nine, nine.nr.int, "0.92", 1)
eleven.rna <- assign_RNA_types_y(eleven, eleven.nr.int, "2.66", 2)
one.rna <- assign_RNA_types_y(one, one.nr.int, "5.74", 3)
three.rna <- assign_RNA_types_y(three, three.nr.int, "8.61", 4)

odRNAy <- rbind(nine.rna,eleven.rna,one.rna,three.rna)

#--bacteria data processing-----------------------------------------------------------------------------

#function to create Seurat Object for downstream analysis from counts, features, barcodes, and other metadata
create_seurat_object_from_seq_files <- function(data_dir, sublibrary, ribosome_removal, threshold){
  
  setwd(list.dirs(data_dir))
  
  data <- readMM('UniqueAndMult-Uniform.mtx')
  genes <- read.table('features.tsv')
  genes <- genes$V1
  genes <- convert_gene_names_bacteria(genes)
  barcodes <- read.table('barcodes.tsv')
  barcodes <- barcodes$V1
  
  rownames(data) <- genes
  colnames(data) <- barcodes
  
  #use filtered cells

  filteredData <- round_one_bc_collapse(data, threshold)
  
  #deal with ribosomal reads
  
  if (ribosome_removal == 'rRNA') {
    ribosomes <- read.csv('C:/Users/lbrettne/Desktop/bacteria.ribosomes.csv')
    colnames(ribosomes) <- c('gene','type')
    rRNA <- ribosomes[which(ribosomes$type == 'rRNA'),]
    mRNA <- filteredData[which(!(rownames(filteredData) %in% rRNA$gene)),]
    
    #convert to Seurat object
    SO <- CreateSeuratObject(counts = mRNA)
  } else {
    #convert to Seurat object
    SO <- CreateSeuratObject(counts = filteredData, min.cells = 15, min.features = 5)
  }
  
  #add metadata
  wells <- assign_cell_wells(colnames(SO@assays$RNA))
  celldata <- wells$well
  celldata <- as.data.frame(celldata)
  row.names(celldata) <- wells$barcode
  
  SO$well <- celldata$celldata
  SO$sublibrary <- sublibrary
  
  return(SO)
}

#B. subtilis replicate 1 
GEO661 <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO661 Solo.out/bacillus_only/GeneFull/raw/','M15','none',0.85)
GEO661.nr <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO661 Solo.out/bacillus_only/GeneFull/raw/','M15','rRNA',0.85)

#B. subtilis replicate 2
GEO660 <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO660 Solo.out/Bacillus Only/GeneFull/raw/','M14','none',0.875)
GEO660.nr <- create_seurat_object_from_seq_files('C:/Users/lbrettne/Desktop/GEO660 Solo.out/Bacillus Only/GeneFull/raw/','M14','rRNA',0.875)

#functions adds cell metadata based on which well barcodes are detected
assign_metadata_bacteria <- function(raw.genefull) {

  s1 <- merge(subset(raw.genefull, subset = well == 'A1'), y = c(subset(raw.genefull, subset = well == 'A2'),
                                                                 subset(raw.genefull, subset = well == 'A3'),
                                                                 subset(raw.genefull, subset = well == 'A4'),
                                                                 subset(raw.genefull, subset = well == 'A5'),
                                                                 subset(raw.genefull, subset = well == 'A6')))
  s1$cond <- 'OD0.5'
  
  s2 <- merge(subset(raw.genefull, subset = well == 'A7'), y = c(subset(raw.genefull, subset = well == 'A8'),
                                                                 subset(raw.genefull, subset = well == 'A9'),
                                                                 subset(raw.genefull, subset = well == 'A10'),
                                                                 subset(raw.genefull, subset = well == 'A11'),
                                                                 subset(raw.genefull, subset = well == 'A12')))
  s2$cond <- 'OD1.0'
  
  s3 <- merge(subset(raw.genefull, subset = well == 'B1'), y = c(subset(raw.genefull, subset = well == 'B2'),
                                                                 subset(raw.genefull, subset = well == 'B3'),
                                                                 subset(raw.genefull, subset = well == 'B4'),
                                                                 subset(raw.genefull, subset = well == 'B5'),
                                                                 subset(raw.genefull, subset = well == 'B6')))
  s3$cond <- 'OD1.3'
  
  s4 <- merge(subset(raw.genefull, subset = well == 'B7'), y = c(subset(raw.genefull, subset = well == 'B8'),
                                                                 subset(raw.genefull, subset = well == 'B9'),
                                                                 subset(raw.genefull, subset = well == 'B10'),
                                                                 subset(raw.genefull, subset = well == 'B11'),
                                                                 subset(raw.genefull, subset = well == 'B12')))
  s4$cond <- 'OD1.6'
  
  s5 <- merge(subset(raw.genefull, subset = well == 'C1'), y = c(subset(raw.genefull, subset = well == 'C2'),
                                                                 subset(raw.genefull, subset = well == 'C3'),
                                                                 subset(raw.genefull, subset = well == 'C4'),
                                                                 subset(raw.genefull, subset = well == 'C5'),
                                                                 subset(raw.genefull, subset = well == 'C6')))
  s5$cond <- 'OD2.8'
  
  s6 <- merge(subset(raw.genefull, subset = well == 'C7'), y = c(subset(raw.genefull, subset = well == 'C8'),
                                                                 subset(raw.genefull, subset = well == 'C9'),
                                                                 subset(raw.genefull, subset = well == 'C10'),
                                                                 subset(raw.genefull, subset = well == 'C11'),
                                                                 subset(raw.genefull, subset = well == 'C12')))
  s6$cond <- 'OD3.6'
  
  s7 <- merge(subset(raw.genefull, subset = well == 'D1'), y = c(subset(raw.genefull, subset = well == 'D2'),
                                                                 subset(raw.genefull, subset = well == 'D3'),
                                                                 subset(raw.genefull, subset = well == 'D4'),
                                                                 subset(raw.genefull, subset = well == 'D5'),
                                                                 subset(raw.genefull, subset = well == 'D6')))
  s7$cond <- 'OD5.3'
  
  s8 <- merge(subset(raw.genefull, subset = well == 'D7'), y = c(subset(raw.genefull, subset = well == 'D8'),
                                                                 subset(raw.genefull, subset = well == 'D9'),
                                                                 subset(raw.genefull, subset = well == 'D10'),
                                                                 subset(raw.genefull, subset = well == 'D11'),
                                                                 subset(raw.genefull, subset = well == 'D12')))
  s8$cond <- 'OD6.0'
  
  growcurv <- merge(s1,y = c(s2,s3,s4,s5,s6,s7,s8))
  
  return(growcurv)
}
assign_metadata_bacteria_M14 <- function(raw.genefull) {
  
  s1 <- merge(subset(raw.genefull, subset = well == 'A1'), y = c(subset(raw.genefull, subset = well == 'A2'),
                                                                 subset(raw.genefull, subset = well == 'A3'),
                                                                 subset(raw.genefull, subset = well == 'A4'),
                                                                 subset(raw.genefull, subset = well == 'A5'),
                                                                 subset(raw.genefull, subset = well == 'A6'),
                                                                 subset(raw.genefull, subset = well == 'A7'),
                                                                 subset(raw.genefull, subset = well == 'A8')))
  s1$cond <- 'OD0.5'
  
  s2 <- merge(subset(raw.genefull, subset = well == 'A9'), y = c(subset(raw.genefull, subset = well == 'A10'),
                                                                 subset(raw.genefull, subset = well == 'A11'),
                                                                 subset(raw.genefull, subset = well == 'A12'),
                                                                 subset(raw.genefull, subset = well == 'B1'),
                                                                 subset(raw.genefull, subset = well == 'B2'),
                                                                 subset(raw.genefull, subset = well == 'B3'),
                                                                 subset(raw.genefull, subset = well == 'B4')))
  s2$cond <- 'OD1.0'
  
  s3 <- merge(subset(raw.genefull, subset = well == 'B5'), y = c(subset(raw.genefull, subset = well == 'B6'),
                                                                 subset(raw.genefull, subset = well == 'B7'),
                                                                 subset(raw.genefull, subset = well == 'B8'),
                                                                 subset(raw.genefull, subset = well == 'B9'),
                                                                 subset(raw.genefull, subset = well == 'B10'),
                                                                 subset(raw.genefull, subset = well == 'B10'),
                                                                 subset(raw.genefull, subset = well == 'B12')))
  s3$cond <- 'OD1.7'
  
  s4 <- merge(subset(raw.genefull, subset = well == 'C1'), y = c(subset(raw.genefull, subset = well == 'C2'),
                                                                 subset(raw.genefull, subset = well == 'C3'),
                                                                 subset(raw.genefull, subset = well == 'C4'),
                                                                 subset(raw.genefull, subset = well == 'C5'),
                                                                 subset(raw.genefull, subset = well == 'C6'),
                                                                 subset(raw.genefull, subset = well == 'C7'),
                                                                 subset(raw.genefull, subset = well == 'C8')))
  s4$cond <- 'OD2.0'
  
  s5 <- merge(subset(raw.genefull, subset = well == 'C9'), y = c(subset(raw.genefull, subset = well == 'C10'),
                                                                 subset(raw.genefull, subset = well == 'C11'),
                                                                 subset(raw.genefull, subset = well == 'C12'),
                                                                 subset(raw.genefull, subset = well == 'D1'),
                                                                 subset(raw.genefull, subset = well == 'D2'),
                                                                 subset(raw.genefull, subset = well == 'D3'),
                                                                 subset(raw.genefull, subset = well == 'D4')))
  s5$cond <- 'OD2.8'
  
  s6 <- merge(subset(raw.genefull, subset = well == 'D5'), y = c(subset(raw.genefull, subset = well == 'D6'),
                                                                 subset(raw.genefull, subset = well == 'D7'),
                                                                 subset(raw.genefull, subset = well == 'D8'),
                                                                 subset(raw.genefull, subset = well == 'D9'),
                                                                 subset(raw.genefull, subset = well == 'D10'),
                                                                 subset(raw.genefull, subset = well == 'D11'),
                                                                 subset(raw.genefull, subset = well == 'D12')))
  s6$cond <- 'OD3.2'
  
  growcurv <- merge(s1,y = c(s2,s3,s4,s5,s6))
  return(growcurv)
}

growcurvb <- assign_metadata_bacteria(GEO661)
growcurvb.nr <- assign_metadata_bacteria(GEO661.nr)

#subset data into individual timepoints
OD0.5 <- subset(growcurvb, subset = cond == "OD0.5")
OD1.0 <- subset(growcurvb, subset = cond == "OD1.0")
OD1.3 <- subset(growcurvb, subset = cond == "OD1.3")
OD1.6 <- subset(growcurvb, subset = cond == "OD1.6")
OD2.8 <- subset(growcurvb, subset = cond == "OD2.8")
OD3.6 <- subset(growcurvb, subset = cond == "OD3.6")
OD5.3 <- subset(growcurvb, subset = cond == "OD5.3")
OD6.0 <- subset(growcurvb, subset = cond == "OD6.0")

#subset data minus rRNA into individual timepoints
OD0.5.nr <- subset(growcurvb.nr, subset = cond == "OD0.5")
OD1.0.nr <- subset(growcurvb.nr, subset = cond == "OD1.0")
OD1.3.nr <- subset(growcurvb.nr, subset = cond == "OD1.3")
OD1.6.nr <- subset(growcurvb.nr, subset = cond == "OD1.6")
OD2.8.nr <- subset(growcurvb.nr, subset = cond == "OD2.8")
OD3.6.nr <- subset(growcurvb.nr, subset = cond == "OD3.6")
OD5.3.nr <- subset(growcurvb.nr, subset = cond == "OD5.3")
OD6.0.nr <- subset(growcurvb.nr, subset = cond == "OD6.0")

#normalize and scale mRNA datasets
scale_normalize_seurat_data <- function(seuratobject, nvariableFeatures, npcs, features){
  
  integrated <- NormalizeData(seuratobject, verbose = FALSE)
  integrated <- FindVariableFeatures(integrated, selection.method = "vst", 
                                     nfeatures = nvariableFeatures, verbose = FALSE)
  
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE, features = features)
  #integrated <- JackStraw(integrated)
  #integrated <- ScoreJackStraw(integrated, dims = 1:npcs)
  #p <- plot_grid(ElbowPlot(integrated), JackStrawPlot(integrated, dims = 1:npcs))
  p <- ElbowPlot(integrated)
  print(p)
  return(integrated)
}

OD0.5.nr.int <- scale_normalize_seurat_data(OD0.5.nr, 3000, 20, NULL)
OD1.0.nr.int <- scale_normalize_seurat_data(OD1.0.nr, 3000, 20, NULL)
OD1.3.nr.int <- scale_normalize_seurat_data(OD1.3.nr, 3000, 20, NULL)
OD1.6.nr.int <- scale_normalize_seurat_data(OD1.6.nr, 3000, 20, NULL)
OD2.8.nr.int <- scale_normalize_seurat_data(OD2.8.nr, 3000, 20, NULL)
OD3.6.nr.int <- scale_normalize_seurat_data(OD3.6.nr, 3000, 20, NULL)
OD5.3.nr.int <- scale_normalize_seurat_data(OD5.3.nr, 3000, 20, NULL)
OD6.0.nr.int <- scale_normalize_seurat_data(OD6.0.nr, 3000, 20, NULL)

#calculate counts and statistics for RNA types, including rRNA, mRNA, etc.
assign_RNA_types_b <- function(SO, SO2, group_name, timepoint) {
  
  ribosomes <- read.csv('C:/Users/lbrettne/Desktop/bacteria.ribosomes.csv')
  colnames(ribosomes) = c("name","type")
  rRNA.genes <- ribosomes[which(ribosomes$type == 'rRNA'),]
  tRNA.genes <- ribosomes[which(ribosomes$type == 'tRNA'),]
  rProt.genes <- ribosomes[which(ribosomes$type == 'rProtein'),]

  class <- read.csv('C:/Users/lbrettne/Desktop/bacteria_RP_iESR_RiBi.csv')
  RP.genes <- class[which(class$class == 'RP'),]
  iESR.genes <- class[which(class$class == 'iESR'),]
  RiBi.genes <- class[which(class$class == 'RiBi'),]
  
  RP  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RP.genes$gene),])
  rpcv <- sd(RP[which(RP>0)])/mean(RP[which(RP>0)])
  iESR  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% iESR.genes$gene),])
  RiBi  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RiBi.genes$gene),])
  totalmRNA <- colSums(SO2@assays$RNA)
  
  SO2$rsect <- 0
  SO2$psect <- 0
  
  for (i in 1:dim(SO2)[2]) {
    
    SO2$rsect[i] <- (RP[i]+RiBi[i])/(RP[i]+RiBi[i]+iESR[i])
    SO2$psect[i] <- (iESR[i])/(RP[i]+RiBi[i]+iESR[i])
    
    
  }
  
  dt <- c(0,1.3509309,1.2972221,1.2699011,1.2351942,1.2072780,1.1572952,1.1383160,1.0740787,0.9988778,0.9714859,0.8519719,0.7547571,0.6883084,0.3290760,0.1606523)
  
  rRNA  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rRNA.genes$name),])
  tRNA  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% tRNA.genes$name),])
  mRNA  <- SO@assays$RNA[which(!(rownames(SO@assays$RNA) %in% rRNA.genes$name)),]
  #mRNA  <- SO[which(!(rownames(mRNA) %in% tRNA.genes$name)),]
  mRNA  <- colSums(mRNA)
  rProt <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rProt.genes$name),])
  rProt.ratio <- rProt/mRNA
  totalRNA <- colSums(SO@assays$RNA)
  rRNA.ratio <- rRNA/totalRNA
  mRNA.ratio <- mRNA/totalRNA
  mRNA.rRNA <- mRNA/rRNA
  bulkrr <- sum(rRNA)/sum(totalRNA)

  
  bulkrp <- log10(sum(rProt)/sum(totalRNA)*10^6)
  rcv <- sd(rRNA)/mean(rRNA)
  rrcv <- sd(rRNA.ratio)/mean(rRNA.ratio)
  #rpcv <- sd(rProt[which(rProt>0)])/mean(rProt[which(rProt>0)])
  rr.med.dev.med <- median(abs(rRNA - median(rRNA)))
  
  temp <- data.frame(group_name = group_name, 
                     rRNA = rRNA, 
                     tRNA = tRNA, 
                     mRNA = mRNA, 
                     rProt = rProt,
                     totalRNA = totalRNA,
                     rProt.ratio = rProt.ratio,
                     rRNA.ratio = rRNA.ratio,
                     mRNA.ratio = mRNA.ratio,
                     mRNA.rRNA = mRNA.rRNA,
                     bulkrr = bulkrr,
                     #GpH = 1/(dt[timepoint]),
                     GpH = dt[timepoint],
                     bulkrp = bulkrp,
                     rcv = rcv,
                     rpcv = rpcv,
                     rrcv = rrcv,
                     rsd = sd(rRNA),
                     rrsd = sd(rRNA.ratio),
                     rpsd = sd(rProt),
                     rprsd = sd(rProt.ratio),
                     rrm = mean(rRNA),
                     rr.med.dev.med = rr.med.dev.med,
                     rsect = SO2$rsect,
                     psect = SO2$psect)
  
  return(temp)
}

OD0.5.RNA <- assign_RNA_types_b(OD0.5, OD0.5.nr.int, "0.5", 4)
OD1.0.RNA <- assign_RNA_types_b(OD1.0, OD1.0.nr.int, "1.0", 5)
OD1.3.RNA <- assign_RNA_types_b(OD1.3, OD1.3.nr.int, "1.3", 8)
OD1.6.RNA <- assign_RNA_types_b(OD1.6, OD1.6.nr.int, "1.6", 9)
OD2.8.RNA <- assign_RNA_types_b(OD2.8, OD2.8.nr.int, "2.8", 12)
OD3.6.RNA <- assign_RNA_types_b(OD3.6, OD3.6.nr.int, "3.6", 14)
OD5.3.RNA <- assign_RNA_types_b(OD5.3, OD5.3.nr.int, "5.3", 15)
OD6.0.RNA <- assign_RNA_types_b(OD6.0, OD6.0.nr.int, "6.0", 16)

odRNAb <- rbind(OD0.5.RNA,OD1.0.RNA,OD1.3.RNA,OD1.6.RNA,OD2.8.RNA,OD3.6.RNA,OD5.3.RNA,OD6.0.RNA)

#repeat for replicate 2
growcurvb.m14 <- assign_metadata_bacteria_M14(GEO660)
growcurvb.nr.m14 <- assign_metadata_bacteria_M14(GEO660.nr)

OD0.5.m14 <- subset(growcurvb.m14, subset = cond == "OD0.5")
OD1.0.m14 <- subset(growcurvb.m14, subset = cond == "OD1.0")
OD1.7.m14 <- subset(growcurvb.m14, subset = cond == "OD1.7")
OD2.0.m14 <- subset(growcurvb.m14, subset = cond == "OD2.0")
OD2.8.m14 <- subset(growcurvb.m14, subset = cond == "OD2.8")
OD3.2.m14 <- subset(growcurvb.m14, subset = cond == "OD3.2")

OD0.5.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD0.5")
OD1.0.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD1.0")
OD1.7.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD1.7")
OD2.0.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD2.0")
OD2.8.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD2.8")
OD3.2.nr.m14 <- subset(growcurvb.nr.m14, subset = cond == "OD3.2")

OD0.5.nr.m14.int <- scale_normalize_seurat_data(OD0.5.nr.m14, 3000, 20, NULL)
OD1.0.nr.m14.int <- scale_normalize_seurat_data(OD1.0.nr.m14, 3000, 20, NULL)
OD1.7.nr.m14.int <- scale_normalize_seurat_data(OD1.7.nr.m14, 3000, 20, NULL)
OD2.0.nr.m14.int <- scale_normalize_seurat_data(OD2.0.nr.m14, 3000, 20, NULL)
OD2.8.nr.m14.int <- scale_normalize_seurat_data(OD2.8.nr.m14, 3000, 20, NULL)
OD3.2.nr.m14.int <- scale_normalize_seurat_data(OD3.2.nr.m14, 3000, 20, NULL)

m14OD0.5.RNA <- assign_RNA_types_b(OD0.5.m14, OD0.5.nr.m14.int, "0.5", 4)
m14OD1.0.RNA <- assign_RNA_types_b(OD1.0.m14, OD1.0.nr.m14.int, "1.0", 5)
m14OD1.7.RNA <- assign_RNA_types_b(OD1.7.m14, OD1.7.nr.m14.int, "1.7", 10)
m14OD2.0.RNA <- assign_RNA_types_b(OD2.0.m14, OD2.0.nr.m14.int, "2.0", 11)
m14OD2.8.RNA <- assign_RNA_types_b(OD2.8.m14, OD2.8.nr.m14.int, "2.8", 12)
m14OD3.2.RNA <- assign_RNA_types_b(OD3.2.m14, OD3.2.nr.m14.int, "3.2", 13)

odRNAb.m14 <- rbind(m14OD0.5.RNA,m14OD1.0.RNA,m14OD1.7.RNA,m14OD2.0.RNA,m14OD2.8.RNA,m14OD3.2.RNA)

#--Figure 1-----------------------------------------------------------------------------

doublingrate <- function(fit, x, t){
  drs <- c()
  for (i in 2:length(fit)){
    dr <- log2(fit[i]/fit[i-1])/(x[i]-x[i-1])
    drs <- c(drs,dr)
  }
  
  grs <- c()
  for (j in 2:length(t)){
    k <- which.min(abs(x - t[j]))
    gr <- drs[k]
    grs <- c(grs, gr)
  }
  return(grs)
}

#yeast
ty <- c(0,     #time in hours
       17.25, #9AM  *sequenced
       19.25, #11AM *
       21.25, #1PM  *
       23.25, #3PM  *
       25.25  #5PM
       )

d1 <- c(13750,8218300,23386700,55516100,87744800,76303100)  #Brettner et al. 2024 (Yeast)
d2 <- c(13750,12568600,38656100,71317800,81984100,88041000) #cells/mL
h1 <- c(13750,9108600,25980300,58423200,87591400,84058800)
h2 <- c(13750,7035600,18339300,44139200,87052400,81699000)

tfity <- rep(ty, 4)
d <- c(d1,d2)
h <- c(h1,h2)
y <- c(d,h)
ymean <- c(mean(d1[2], d2[2], h1[2], h2[2]),
           mean(d1[3], d2[3], h1[3], h2[3]),
           mean(d1[4], d2[4], h1[4], h2[4]),
           mean(d1[5], d2[5], h1[5], h2[5]))

fity <- nls(y ~ (L*x0*exp(k*tfity))/(L+x0*(exp(k*tfity)-1)), start = list(L = 10e08, k = 0.8, x0 = 20))
Ly  <- summary(fity)$parameters[1,1]
ky  <- summary(fity)$parameters[2,1]
x0y <- summary(fity)$parameters[3,1]
xy <- seq(from = 0, to = 26, by = 0.1)
fity <- (Ly*x0y*exp(ky*xy))/(Ly+x0y*(exp(ky*xy)-1))

GpHy <- doublingrate(fity,xy,xy)
GpHy <- data.frame(t = xy[-1], gph = GpHy)
GpHy_data <- doublingrate(fity,xy,ty)
GpHy_data <- data.frame(t = ty[c(-1,-6)], gph = GpHy_data[-5])

cell_density <- data.frame(t = xy, cpml = fity)
cell_dens_data <- data.frame(t = ty[c(-1,-6)], cpml = ymean )

y1.1.1 <- ggplot(GpHy, aes(x = t, y = gph)) +
  geom_line() +
  geom_point(data = GpHy_data, aes(x = t, y = gph), size = 4, fill = "#FFEBC0", shape = 21) +
  ylab("growth rate (generations/hour)") +
  xlab("time (hours)") +
  scale_y_continuous(position = "right") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
y1.1.1

y1.1.2 <- ggplot(cell_density, aes(x = t, y = cpml/10^7)) +
  geom_line(color = "gray") +
  geom_point(data = cell_dens_data, aes(x = t, y = cpml/10^7), size = 4, fill = "#FFEBC0", shape = 21, color = "gray") +
  ylab("growth rate (generations/hour)") +
  xlab("time (hours)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
y1.1.2

y1.2 <- ggplot(odRNAy, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
y1.2

meanrRNAy <- c(exp(mean(log(odRNAy$rRNA[which(odRNAy$group_name == "0.92")]))),
              exp(mean(log(odRNAy$rRNA[which(odRNAy$group_name == "2.66")]))),
              exp(mean(log(odRNAy$rRNA[which(odRNAy$group_name == "5.74")]))),
              exp(mean(log(odRNAy$rRNA[which(odRNAy$group_name == "8.61")]))))
GpHy <- GpHy_data$gph
summary(lm(meanrRNAy~GpHy))

y1.3 <- ggplot(odRNAy, aes(x = GpH, y = log10(rRNA), fill = factor(GpH))) +
  #geom_jitter(aes(color = ploidy), alpha = 0.5, size = 0.05) +
  geom_jitter(color = "#9b870c", alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  #stat_summary(aes(color = ploidy),fun = mean, geom = "point", size = 4, shape = 21) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  #scale_color_manual(values = c("#9b870c","#643B9F")) +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent', color = NA),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()
y1.3

summary(lm(log10(odRNAy$rRNA)~odRNAy$GpH))
summary(lm(odRNAy$rRNA~odRNAy$GpH))

#bacteria
hour  <- c(0.00,1.10,2.13,2.45,2.75,2.90,3.18,3.28,3.60,3.93,4.0,4.42,4.67,4.95,6.12,6.97) #Kuchina and Brettner et al. 2021 (Science)
OD600 <- c(0.02,0.04,0.27,0.49,0.94,0.93,1.16,1.34,1.58,1.76,2.0,2.87,3.20,3.56,5.31,6.34)

fitb <- nls(OD600 ~ L/(1+exp(-k*(hour - x0))), start = list(L = 6, k = 0.5, x0 = 3))
Lb  <- summary(fitb)$parameters[1,1]
kb  <- summary(fitb)$parameters[2,1]
x0b <- summary(fitb)$parameters[3,1]
xb <- seq(from = 0, to = 8, by = 0.1)
fitb <- Lb/(1+exp(-kb*(xb - x0b)))

GpHb <- doublingrate(fitb,xb,xb)
GpHb <- data.frame(t = xb[-1], gph = GpHb)
GpHb_data <- doublingrate(fitb,xb,hour)
GpHb_data <- data.frame(t = hour[c(4,5,8,9,12,13,14,15)], gph = GpHb_data[c(3,4,7,8,11,12,13,14)])

cell_density <- data.frame(t = xb, cpml = fitb)
cell_dens_data <- data.frame(t = hour[c(4,5,8,9,12,13,14,15)], cpml = OD600[c(4,5,8,9,12,13,14,15)] )

b1.1.1 <- ggplot(GpHb, aes(x = t, y = gph)) +
  geom_line() +
  geom_point(data = GpHb_data, aes(x = t, y = gph), size = 4, fill = "#8CA2CB", shape = 21) +
  ylab("growth rate (generations/hour)") +
  xlab("time (hours)") +
  scale_y_continuous(position = "right") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b1.1.1

b1.1.2 <- ggplot(cell_density, aes(x = t, y = cpml)) +
  geom_line(color = "gray") +
  geom_point(data = cell_dens_data, aes(x = t, y = cpml), size = 4, fill = "#8CA2CB", shape = 21, color = "gray") +
  ylab("growth rate (generations/hour)") +
  xlab("time (hours)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b1.1.2

b1.2 <- ggplot(odRNAb, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b1.2

x <- c(mean(OD0.5.RNA$GpH), 
       mean(OD1.0.RNA$GpH), 
       mean(OD1.3.RNA$GpH), 
       mean(OD1.6.RNA$GpH),
       mean(OD2.8.RNA$GpH),
       mean(OD3.6.RNA$GpH),
       mean(OD5.3.RNA$GpH),
       mean(OD6.0.RNA$GpH))

y <- c(mean(OD0.5.RNA$rRNA), 
       mean(OD1.0.RNA$rRNA), 
       mean(OD1.3.RNA$rRNA), 
       mean(OD1.6.RNA$rRNA),
       mean(OD2.8.RNA$rRNA),
       mean(OD3.6.RNA$rRNA),
       mean(OD5.3.RNA$rRNA),
       mean(OD6.0.RNA$rRNA))

summary(lm(y~x))

b1.3 <- ggplot(odRNAb, aes(x = GpH, y = log10(rRNA), fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(2,3,4,5),labels = c("100","1,000", "10,000", "100,000")) +
  NoLegend()
b1.3

summary(lm(log10(odRNAb$rRNA)~odRNAb$GpH))

#--Figure 2-----------------------------------------------------------------------------

y2.1 <- ggplot(odRNAy, aes(x = GpH, y = rcv, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("ribosome variation CV(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend() +
  xlab("population GR") +
  ylab("rRNA counts/cell CV")
y2.1

y2.2 <- ggplot(odRNAy, aes(x = GpH, y = rpcv, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("ribosome variation CV(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
  ) +
  NoLegend() +
  xlab("population GR") +
  ylab("rProtein mRNA CV")
y2.2

b2.1 <- ggplot(odRNAb, aes(x = GpH, y = rcv, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  #stat_summary(data = odRNAb.m14, aes(x = GpH, y = rcv, fill = factor(GpH)), fun = mean, geom = "point", size = 4, fill = "#2c4470", shape = 21) +
  ylab("ribosome variation CV(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend() +
  xlab("population GR") +
  ylab("rRNA counts/cell CV")
b2.1

b2.2 <- ggplot(odRNAb, aes(x = GpH, y = rpcv, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  #stat_summary(data = odRNAb.m14, aes(x = GpH, y = rcv, fill = factor(GpH)), fun = mean, geom = "point", size = 4, fill = "#2c4470", shape = 21) +
  ylab("ribosome variation CV(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
  ) +
  NoLegend() +
  ylab("rProtein mRNA CV") +
  xlab("population GR")
b2.2

scale_normalize_seurat_data <- function(seuratobject, nvariableFeatures, npcs, features){
  
  integrated <- NormalizeData(seuratobject, verbose = FALSE)
  integrated <- FindVariableFeatures(integrated, selection.method = "vst", 
                                     nfeatures = nvariableFeatures, verbose = FALSE)
  
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE, features = features)
  #integrated <- JackStraw(integrated)
  #integrated <- ScoreJackStraw(integrated, dims = 1:npcs)
  #p <- plot_grid(ElbowPlot(integrated), JackStrawPlot(integrated, dims = 1:npcs))
  p <- ElbowPlot(integrated)
  print(p)
  return(integrated)
}

nine.nr.int<- scale_normalize_seurat_data(nine.nr, 3000, 20, NULL)
eleven.nr.int <- scale_normalize_seurat_data(eleven.nr, 3000, 20, NULL)
one.nr.int <- scale_normalize_seurat_data(one.nr, 3000, 20, NULL)
three.nr.int <- scale_normalize_seurat_data(three.nr, 3000, 20, NULL)

nine.int<- scale_normalize_seurat_data(nine, 3000, 20, NULL)
eleven.int <- scale_normalize_seurat_data(eleven, 3000, 20, NULL)
one.int <- scale_normalize_seurat_data(one, 3000, 20, NULL)
three.int <- scale_normalize_seurat_data(three, 3000, 20, NULL)

rRNAnBrauer <- function(SO, SO1, SO2, group_name) {
  
  brauer_gr_slopes <- read.csv('C:/Users/lbrettne/Desktop/brauer_gr_slopes_nozeroes.csv')
  brauer_gr_slopes$ORF <- convert_gene_names_yeast(brauer_gr_slopes$Name, FALSE)
  
  ribosomes <- read.csv('C:/Users/lbrettne/Desktop/Solo.out/ribosomes.csv')
  names(ribosomes)[1] = "name1"
  rRNA.genes <- ribosomes[which(ribosomes$type == 'rRNA'),]
  colnames(rRNA.genes) <- c("name1", "name2", "type", "location")
  tRNA.genes <- ribosomes[which(ribosomes$type == 'tRNA'),]
  colnames(tRNA.genes) <- c("name1", "name2", "type", "location")
  rProt.genes <- ribosomes[which(ribosomes$type == 'rProtein'),]
  colnames(rProt.genes) <- c("name1", "name2", "type", "location")
  
  class <- read.csv('C:/Users/lbrettne/Desktop/RiBi_iESR_CC_genes.csv')
  class$Gene <- convert_gene_names_yeast(class$Gene, FALSE)
  RP.genes <- class[which(class$Group == 'RP'),]
  iESR.genes <- class[which(class$Group == 'iESR'),]
  RiBi.genes <- class[which(class$Group == 'RiBi'),]

  RP  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RP.genes$Gene),])
  iESR  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% iESR.genes$Gene),])
  RiBi  <- colSums(SO2@assays$RNA[which(rownames(SO2@assays$RNA) %in% RiBi.genes$Gene),])
  totalmRNA <- colSums(SO2@assays$RNA)
  
  temp.bgr <- brauer_gr_slopes[which(brauer_gr_slopes$ORF %in% row.names(SO2@assays$RNA)),]
  brauer.down <- temp.bgr[which(temp.bgr$GRR == "down"),]
  brauer.up <- temp.bgr[which(temp.bgr$GRR == "up"),]
  
  temp1.exp <- SO1@assays$RNA[which(row.names(SO1@assays$RNA) %in% brauer_gr_slopes$ORF),]
  temp2.exp <- SO2@assays$RNA[which(row.names(SO2@assays$RNA) %in% brauer_gr_slopes$ORF),]
  temp1.exp <- temp1.exp[temp.bgr$ORF,]
  temp2.exp <- temp2.exp[temp.bgr$ORF,]
  
  bulk = (1/sum(as.matrix(SO2@assays$RNA[])))*(sum(rowSums(temp2.exp)/temp.bgr$Slope))
  dR = sum((log2(rowSums(temp2.exp)/sum(SO2@assays$RNA[])+1) - log2(rowSums(temp1.exp)/sum(SO1@assays$RNA[])+1))/temp.bgr$Slope)
  SO2$gr <- 0
  SO2$up.pie <- 0
  SO2$down.pie <- 0
  
  up_down <- c()
  
  for (i in 1:dim(temp2.exp)[2]) {
    gr <- (1/(sum(as.matrix(SO2@assays$RNA[,i]))))*(sum(temp2.exp[,i]/temp.bgr$Slope))
    #gr <- (sum(temp.exp[,i]/temp.bgr$Slope))
    SO2$gr[i] <- gr
    
    up.pie <- sum(temp2.exp[which(rownames(temp2.exp) %in% brauer.up$ORF),i])/sum(temp2.exp[,i])
    down.pie <- sum(temp2.exp[which(rownames(temp2.exp) %in% brauer.down$ORF),i])/sum(temp2.exp[,i])
    
    SO2$up.pie[i] <- up.pie
    SO2$down.pie[i] <- down.pie

    if(up.pie >= 0.75) {
      up_down[i] = "UP"
    } else if(down.pie >= 0.75) {
      up_down[i] = "DOWN"
    } else {
      up_down[i] = "NONE"
    }
  }
  
  SO2$rsect <- 0
  SO2$psect <- 0
  
  for (i in 1:dim(SO2)[2]) {

    SO2$rsect[i] <- (RP[i]+RiBi[i])/(RP[i]+RiBi[i]+iESR[i])
    SO2$psect[i] <- (iESR[i])/(RP[i]+RiBi[i]+iESR[i])
    

  }
  
  rRNA  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rRNA.genes$name1 | rownames(SO@assays$RNA) %in% rRNA.genes$name2),])
  totalRNA <- colSums(SO@assays$RNA)
  rRNA.ratio <- rRNA/totalRNA
  bulkrr <- log10(sum(rRNA)/sum(totalRNA)*10^6)
  
  Gxt <- GpHy_data$gph
  
  temp <- data.frame(group_name = group_name, 
                     rRNA = rRNA, 
                     totalRNA = totalRNA,
                     rRNA.ratio = rRNA.ratio,
                     bulkrr = bulkrr,
                     gr = SO2$gr,
                     #grcv = sd(SO1$gr)/mean(SO1$gr),
                     GpH = Gxt[group_name],
                     bulk = bulk,
                     dR = dR,
                     up_down = up_down,
                     up = length(which(up_down == "UP"))/length(up_down),
                     down = length(which(up_down == "DOWN"))/length(up_down),
                     up.pie = SO2$up.pie,
                     down.pie = SO2$down.pie,
                     rsect = SO2$rsect,
                     psect = SO2$psect)
  
  return(temp)
} #calculate Brauer estimated growth rates

nine_nine.br <- rRNAnBrauer(nine, nine.nr.int, nine.nr.int, 1)
nine_eleven.br <- rRNAnBrauer(eleven, nine.nr.int, eleven.nr.int, 2)
nine_one.br <- rRNAnBrauer(one, nine.nr.int, one.nr.int, 3)
nine_three.br <- rRNAnBrauer(three, nine.nr.int, three.nr.int, 4)

baseline_nine <- rbind(nine_nine.br, nine_eleven.br, nine_one.br, nine_three.br) 
baseline_nine$dR <- (baseline_nine$dR + abs(min(baseline_nine$dR)))/abs(min(baseline_nine$dR))
baseline_nine$gr_norm <- (baseline_nine$gr + abs(min(baseline_nine$gr)))/abs(min(baseline_nine$gr))

#nine.br <- rRNAnBrauer(nine, nine.nr.int, 1)
#eleven.br <- rRNAnBrauer(eleven, eleven.nr.int, 2)
#one.br <- rRNAnBrauer(one, one.nr.int, 3)
#three.br <- rRNAnBrauer(three, three.nr.int, 4)

nine.br$grcv <- sd(nine.br$gr + abs(min(nine.br$gr)))/mean(nine.br$gr + abs(min(three.br$gr)))
eleven.br$grcv <- sd(eleven.br$gr + abs(min(eleven.br$gr)))/mean(eleven.br$gr + abs(min(three.br$gr)))
one.br$grcv <- sd(one.br$gr + abs(min(one.br$gr)))/mean(one.br$gr + abs(min(three.br$gr)))
three.br$grcv <- sd(three.br$gr + abs(min(three.br$gr)))/mean(three.br$gr + abs(min(three.br$gr)))

all.br <- rbind(nine.br, eleven.br, one.br, three.br)
all.br$bulk <- ((all.br$bulk + max(all.br$bulk))/(2*max(all.br$bulk)))*1.056114

brauer_fit <- nls(grcv ~ A*exp(-b*GpH) + D, data = all.br, start = list(A = 1.6, b = 4, D = 1))

x = seq(from = 0.1002274, to = 1.056114, by = 0.01)
brauer_fit_line <- data.frame(x = x, y = 3.774e-01*exp(-8.053*x) + 9.616e-02)

y2.2 <- ggplot(all.br, aes(x = GpH, y = grcv)) +
  geom_line(data = brauer_fit_line, aes(x = x, y = y)) +
  stat_summary(aes(fill = factor(GpH)), fun = mean, geom = "point", size = 4, shape = 21) +
  ylab("inferred single-cell growth rate variation CV(BEGR)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#89BB87","#89BB87","#89BB87","#89BB87")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
y2.2

ggplot(baseline_nine, aes(x = GpH, y = log10(rRNA), fill = factor(GpH))) +
  geom_jitter(aes(color = up_down), alpha = 0.5, size = 0.5) +
  geom_violin(trim = FALSE, alpha = 0) +
  #stat_summary(data = diploids, aes(x = GpH, y = log10(rRNA), fill = ploidy),fun = mean, geom = "point", size = 4, shape = 21) +
  #stat_summary(data = haploids, aes(x = GpH, y = log10(rRNA), fill = ploidy),fun = mean, geom = "point", size = 4, shape = 21) +
  #geom_smooth(data = diploids, aes(x = GpH, y = log10(rRNA), fill = ploidy, color = ploidy), method = "lm", se = FALSE, size = 1.5) +
  #geom_smooth(data = haploids, aes(x = GpH, y = log10(rRNA), fill = ploidy, color = ploidy), method = "lm", se = FALSE, size = 1.5, linetype = "dashed") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_color_manual(values = c("red","gray", "green")) +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0","#9b870c","#643B9F")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)#,
    #axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()

ys2.1 <- ggplot(baseline_nine, aes(x = GpH, y = gr_norm, fill = factor(GpH))) +
  geom_jitter(color = "#99CCB2", alpha = 0.5, size = 0.01) +
  geom_violin(position = "dodge", trim = FALSE, alpha = 0.8) +
  #stat_summary(aes(x = GpH, y = dR), fun = mean, geom = "point", size = 2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("inferred scGR (arbitrary)") +
  xlab("pop GR (generations/hour)") +
  scale_fill_manual(values = c("#89BB87","#89BB87","#89BB87","#89BB87")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
ys2.1

sd(summary(lm(baseline_nine$gr_norm~baseline_nine$GpH))$residuals)

length(which(all.br$gr>0))/length(all.br$gr)
length(which(nine.br$gr>0))/length(nine.br$gr)
length(which(eleven.br$gr>0))/length(eleven.br$gr)
length(which(one.br$gr>0))/length(one.br$gr)
length(which(three.br$gr>0))/length(three.br$gr)

area <- read.table("C:/Users/lbrettne/Desktop/LWOdata/LWOdata/LWOarea.txt", header = TRUE)
data <- read.table("C:/Users/lbrettne/Desktop/LWOdata/LWOdata/LWOdata.txt", header = TRUE)

labhigluc = 1:dim(area)[1]

growthrates <- matrix(nrow = length(labhigluc), ncol = dim(area)[2]-1)

k = 0
for (i in labhigluc){
  temp <- area[i,]
  k = k+1
  
  for (j in 1:(length(temp)-1)){
    if (temp[j] == 0){
      growthrates[k,j] = NA
    } else if (temp[j+1] == 0){
      growthrates[k,j] = NA
    } else {
      growthrates[k,j] = as.numeric(log2(temp[j+1]/temp[j]))
    }
  }
}

grnoise <- data.frame(mean = c(), sd = c(), lower = c(), upper = c(), cv = c(), timepoint = c())

for (i in 1:dim(growthrates)[2]){
  grs <- growthrates[,i]
  grs <- grs[which(grs != "NA")]
  grs <- grs[which(grs != -Inf)]
  grs <- grs[which(grs >= 0)]
  
  temp <- data.frame(mean = mean(grs), sd = sd(grs), lower = mean(grs)-sd(grs), upper = mean(grs)+sd(grs), cv = abs(sd(grs)/mean(grs)), timepoint = i)
  grnoise <- rbind(grnoise, temp)
}

boxplot(growthrates)

plot(grnoise$timepoint, grnoise$mean)
plot(grnoise$timepoint, grnoise$sd)
plot(grnoise$timepoint, grnoise$cv)
plot(grnoise$mean, grnoise$sd)
plot(grnoise$mean, grnoise$cv)

ziv_fit <- nls(cv ~ A*exp(-b*mean) + D, data = grnoise, start = list(A = 1.6, b = 4, D = 1))

x = seq(from = 0.2, to = 0.5, by = 0.01)
ziv_fit_line <- data.frame(x = x, y = 38.79941*exp(-17.21114*x) + 0.49978)

y2.3 <- ggplot(grnoise, aes(x = mean, y = cv)) +
  geom_line(data = ziv_fit_line, aes(x = x, y = y), color = "black") +
  geom_point(size = 4, shape = 21, fill = "#E6E6FA") +
  xlab("mean doubling rate") +
  ylab("CV(single colony doubling rate)") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
y2.3

y2.4 <- ggplot(all.br, aes(x = gr, y = lo10(rRNA))) +
  geom_point(size = 1, colour = "#89BB87") +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("rRNA counts/cell") +
  xlab("inferred single-cell growth rate BEGR") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
y2.4

summary(lm(log10(all.br$rRNA~all).br$gr))



#--Figure 3---------------------------------------------------------------------
scale_normalize_seurat_data <- function(seuratobject, nvariableFeatures, npcs, features){
  
  integrated <- NormalizeData(seuratobject, verbose = FALSE)
  integrated <- FindVariableFeatures(integrated, selection.method = "vst", 
                                     nfeatures = nvariableFeatures, verbose = FALSE)
  
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE, features = features)
  #integrated <- JackStraw(integrated)
  #integrated <- ScoreJackStraw(integrated, dims = 1:npcs)
  #p <- plot_grid(ElbowPlot(integrated), JackStrawPlot(integrated, dims = 1:npcs))
  p <- ElbowPlot(integrated)
  print(p)
  return(integrated)
}
find_neighbors_clusters <- function(integratedSO, dims, resolution){
  
  integrated <- FindNeighbors(integratedSO, dims = 1:dims)
  clustered <- FindClusters(integrated, resolution = resolution)
  clustered <- RunUMAP(clustered, dims = 1:dims)
  
  return(clustered)
}
differentially_expressed_genes <- function(clusteredSO) {
  
  markers <- FindAllMarkers(clusteredSO, only.pos = TRUE)
  
  for (j in 0:(length(levels(markers$cluster))-1))
  {
    cluster.genes <- markers[which(markers$cluster == j),]
    cluster.genes <- cluster.genes[which(cluster.genes$p_val_adj < 0.05),]
    print(head(cluster.genes[order(-cluster.genes$avg_log2FC),], n = 100))
  }
}

nine.nr.clust <- find_neighbors_clusters(nine.nr.int, 10, 0.2)
nine.nr.clust$rRNA <- nine.rna$rRNA

three.nr.clust <- find_neighbors_clusters(three.nr.int, 10, 0.1)
three.nr.clust$rRNA <- three.rna$rRNA
three.nr.clust$seurat_clusters <- as.numeric(three.nr.clust$seurat_clusters)
three.nr.clust$seurat_clusters[which(three.nr.clust$seurat_clusters == 1)] = 3
three.nr.clust$seurat_clusters[which(three.nr.clust$seurat_clusters == 2)] = 0
three.nr.clust$seurat_clusters[which(three.nr.clust$seurat_clusters == 3)] = 1
three.nr.clust$seurat_clusters <- as.factor(three.nr.clust$seurat_clusters)

y3.1 <- DimPlot(nine.nr.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
      NoLegend() + 
      theme(aspect.ratio = 1) +
      scale_color_manual(values = c("#339966","#99CCB2", "#808080"))
y3.1

y3.2 <- DimPlot(three.nr.clust, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("#99CCB2","#808080")) +
  ggtitle("")
y3.2

differentially_expressed_genes(nine.nr.clust)
differentially_expressed_genes(three.nr.clust)

nine.zero <- subset(nine.nr.clust, subset = seurat_clusters == 0)
nine.one <- subset(nine.nr.clust, subset = seurat_clusters == 1)
nine.two <- subset(nine.nr.clust, subset = seurat_clusters == 2)

nine.02 <- rbind(data.frame(rRNA = nine.zero$rRNA, cluster = 0), 
                 data.frame(rRNA = nine.two$rRNA, cluster = 2))

three.zero <- subset(three.nr.clust, subset = seurat_clusters == 0)
three.one <- subset(three.nr.clust, subset = seurat_clusters == 1)

three.01 <- rbind(data.frame(rRNA = three.zero$rRNA, cluster = 0),
                  data.frame(rRNA = three.one$rRNA, cluster = 1))

FeaturePlot(nine.nr.clust, features = "TSL1", reduction = "umap") + theme(aspect.ratio = 1)

length(which(nine.nr@assays$RNA[which(rownames(nine.nr@assays$RNA) == "TSL1"),] > 1))/3043


y3.3 <- ggplot(nine.02, aes(x = as.factor(cluster), y = log10(rRNA), fill = as.factor(cluster))) +
        geom_jitter(alpha = 0.05) +
        #geom_boxplot(notch = TRUE, alpha = 0.8) +
        geom_violin(trim = FALSE, alpha = 0.8) +
        stat_summary(fun = mean, geom = "point", size = 3) +
        scale_fill_manual(values = c("#339966","#808080")) +
        theme_classic() +
        ylab("rRNA counts/cell") +
        xlab("cluster") +
        theme(
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          axis.line = element_line(size = 0.5, color = "black"),
          aspect.ratio = 1,
          text = element_text(size=16)#,
          #axis.title = element_blank()
        ) +
        scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
        NoLegend()
y3.3

y3.4 <- ggplot(three.01, aes(x = as.factor(cluster), y = log10(rRNA), fill = as.factor(cluster))) +
  geom_jitter(alpha = 0.05) +
  #geom_boxplot(notch = TRUE, alpha = 0.8) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_fill_manual(values = c("#99CCB2","#808080")) +
  theme_classic() +
  ylab("rRNA counts/cell") +
  xlab("cluster") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)#,
    #axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()
y3.4

t.test(log10(nine.zero$rRNA),log10(nine.two$rRNA))

differentially_expressed_genes(nine.nr.clust)

#bacteria 
OD1.0.nr.int <- scale_normalize_seurat_data(OD1.0.nr, 3000, 20, NULL)
OD1.0.nr.clust <- find_neighbors_clusters(OD1.0.nr.int, 10, 0.1)
OD1.0.nr.clust$rRNA <- OD1.0.RNA$rRNA

OD1.3.nr.int <- scale_normalize_seurat_data(OD1.3.nr, 3000, 20, NULL)
OD1.3.nr.clust <- find_neighbors_clusters(OD1.3.nr.int, 10, 0.1)
OD1.3.nr.clust$rRNA <- OD1.3.RNA$rRNA

b3.1 <- DimPlot(OD1.0.nr.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("#339966", "#808080"))
b3.1

b3.2 <- DimPlot(OD1.3.nr.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("#339966", "#808080"))
b3.2

OD1.0.zero <- subset(OD1.0.nr.clust, subset = seurat_clusters == 0)
OD1.0.one <- subset(OD1.0.nr.clust, subset = seurat_clusters == 1)

OD1.0.01 <- rbind(data.frame(rRNA = OD1.0.zero$rRNA, cluster = 0), 
                 data.frame(rRNA = OD1.0.one$rRNA, cluster = 1))

OD1.3.zero <- subset(OD1.3.nr.clust, subset = seurat_clusters == 0)
OD1.3.one <- subset(OD1.3.nr.clust, subset = seurat_clusters == 1)

OD1.3.01 <- rbind(data.frame(rRNA = OD1.3.zero$rRNA, cluster = 0), 
                  data.frame(rRNA = OD1.3.one$rRNA, cluster = 1))

b3.2 <- ggplot(OD1.0.01, aes(x = as.factor(cluster), y = log10(rRNA), fill = as.factor(cluster))) +
  geom_jitter(alpha = 0.05) +
  #geom_boxplot(notch = TRUE, alpha = 0.8) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_fill_manual(values = c("#339966", "#808080")) +
  theme_classic() +
  ylab("rRNA counts/cell") +
  xlab("cluster") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)#,
    #axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()
b3.2

b3.4 <- ggplot(OD1.3.01, aes(x = as.factor(cluster), y = log10(rRNA), fill = as.factor(cluster))) +
  geom_jitter(alpha = 0.05) +
  #geom_boxplot(notch = TRUE, alpha = 0.8) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_fill_manual(values = c("#339966", "#808080")) +
  theme_classic() +
  ylab("rRNA counts/cell") +
  xlab("cluster") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)#,
    #axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()
b3.4

t.test(log10(OD1.0.zero$rRNA),log10(OD1.0.one$rRNA))

differentially_expressed_genes(OD1.0.nr.clust)
differentially_expressed_genes(OD1.3.nr.clust)

#--Figure 4---------------------------------------------------------------------

rRNAquarts <- function(odRNA){
  GpH <- levels(as.factor(odRNA$GpH))
  
  mins <- c()
  maxs <- c()
  
  for (i in 1:length(GpH)) {
    rRNA <- odRNA$rRNA[which(odRNA$GpH == GpH[i])]
    
    mins[i] <- min(rRNA)
    maxs[i] <- max(rRNA)
  }
  
  maxmins = max(mins)
  minmaxs = min(maxs)
  
  rRNArange <- odRNA$rRNA[which(odRNA$rRNA >= maxmins & odRNA$rRNA <= minmaxs)]
  names.r1 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[1] & odRNA$rRNA < summary(rRNArange)[2]),]), quartile = 1)
  names.r2 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[2] & odRNA$rRNA < summary(rRNArange)[3]),]), quartile = 2)
  names.r3 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[3] & odRNA$rRNA < summary(rRNArange)[5]),]), quartile = 3)
  names.r4 <- data.frame(names = row.names(odRNA[which(odRNA$rRNA >= summary(rRNArange)[5] & odRNA$rRNA < summary(rRNArange)[6]),]), quartile = 4)
  
  quartile_bcs <- rbind(names.r1, names.r2, names.r3, names.r4)
  
  return(quartile_bcs)
}

#yeast
quartiles <- rRNAquarts(odRNAy)

yQ1 <- growcurvy.nr[,which(colnames(growcurvy.nr@assays$RNA) %in% quartiles$names[which(quartiles$quartile == 1)])]
yQ4 <- growcurvy.nr[,which(colnames(growcurvy.nr@assays$RNA) %in% quartiles$names[which(quartiles$quartile == 4)])]

yQ1.int <- scale_normalize_seurat_data(yQ1, 3000, 20, NULL)
yQ4.int <- scale_normalize_seurat_data(yQ4, 3000, 20, NULL)

yQ1.clustered <- find_neighbors_clusters(yQ1.int, 10, 0.2)
yQ4.clustered <- find_neighbors_clusters(yQ4.int, 10, 0.1)

y4.1 <- DimPlot(yQ4.clustered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#339966", "#99CCB2", "#808080"))
y4.1

FeaturePlot(yQ4.clustered, features = "Q0020", reduction = "umap")

#differentially_expressed_genes(yQ4.clustered)

y4.2 <-DimPlot(yQ4.clustered, reduction = "umap", group.by = "hour", label = FALSE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) + 
  labs(title = NULL)
y4.2

y4.3 <- DimPlot(yQ1.clustered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#339966", "#808080"))
y4.3

#differentially_expressed_genes(yQ1.clustered)

y4.4 <-DimPlot(yQ1.clustered, reduction = "umap", group.by = "hour", label = FALSE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) + 
  labs(title = NULL)
y4.4

#bacteria
quartiles <- rRNAquarts(odRNAb)

bQ1 <- growcurvb.nr[,which(colnames(growcurvb.nr@assays$RNA) %in% quartiles$names[which(quartiles$quartile == 1)])]
bQ4 <- growcurvb.nr[,which(colnames(growcurvb.nr@assays$RNA) %in% quartiles$names[which(quartiles$quartile == 4)])]

bQ1.int <- scale_normalize_seurat_data(bQ1, 3000, 20, NULL)
bQ4.int <- scale_normalize_seurat_data(bQ4, 3000, 20, NULL)

bQ1.clustered <- find_neighbors_clusters(bQ1.int, 10, 0.2)
bQ4.clustered <- find_neighbors_clusters(bQ4.int, 10, 0.1)

b4.1 <- DimPlot(bQ4.clustered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#339966", "#99CCB2", "#808080"))
b4.1

#differentially_expressed_genes(bQ4.clustered)

b4.2 <-DimPlot(bQ4.clustered, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#B76E79","#8A3324","#ED820E","#FFD700","#C7E065","#006666","#90E0EF","#B47EDE")) + 
  labs(title = NULL)
b4.2

b4.3 <- DimPlot(bQ1.clustered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#339966", "#808080","#99CCB2"))
b4.3

#differentially_expressed_genes(bQ1.clustered)

b4.4 <-DimPlot(bQ1.clustered, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1) + 
  scale_color_manual(values = c("#B76E79","#8A3324","#ED820E","#FFD700","#C7E065","#006666","#90E0EF","#B47EDE")) + 
  labs(title = NULL)
b4.4

plot_grid(y4.1,b4.1,y4.3,b4.3)
plot_grid(y4.2,b4.2,y4.4,b4.4)





#--Figure 5---------------------------------------------------------------------

b5.1 <- ggplot(test, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
b5.1

b5.2 <- ggplot(test, aes(x = GpH, y = log10(mRNA), fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(2,3,4,5),labels = c("100","1,000", "10,000", "100,000")) +
  NoLegend()
b5.2

growcurvb.nr.m14.int <- scale_normalize_seurat_data(growcurvb.nr.m14, 3000, 20, NULL)
growcurvb.nr.m14.clust <- find_neighbors_clusters(growcurvb.nr.m14.int, 10, 0.05)
growcurvb.nr.m14.clust$mRNA.rRNA <- odRNAb.m14$mRNA.rRNA

b5.3 <- DimPlot(growcurvb.nr.m14.clust, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, pt.size = 0.01) + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("#B76E79","#8A3324","#ED820E","#FFD700","#C7E065","#006666"))+ NoLegend()
b5.3

DimPlot(growcurvb.nr.m14.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1)

b5.4 <- FeaturePlot(growcurvb.nr.m14.clust, reduction = "umap", features = "mRNA.rRNA", cols = c("black","#8CA2CB", "#89BB87"), pt.size = 0.01)+ 
  theme(aspect.ratio = 1) + NoLegend()
b5.4

b5.3+b5.4

differentially_expressed_genes(growcurvb.nr.m14.clust)

spores <- subset(growcurvb.nr.m14.clust, subset = seurat_clusters == 2)
nsspores <- subset(growcurvb.nr.m14.clust, subset = seurat_clusters != 2)

test0 <- subset(nine.nr.clust, subset = seurat_clusters == 0) 
test1 <- subset(nine.nr.clust, subset = seurat_clusters == 1)

test <- merge(test0, y = test1)
test.int <- scale_normalize_seurat_data(test, 3000, 20, NULL)
test.clust <- find_neighbors_clusters(test.int, 10, 0.2)

differentially_expressed_genes(test.clust)

#--Supplemental Figures---------------------------------------------------------
#--SFs additional means and distributions-------------------------------------------------------------------------

#yeast

nine.int.rna <- assign_RNA_types_y(nine.int,1,1)
eleven.int.rna <- assign_RNA_types_y(eleven.int,2,2)
one.int.rna <- assign_RNA_types_y(one.int,3,3)
three.int.rna <- assign_RNA_types_y(three.int,4,4)

odRNAy.int <- rbind(nine.int.rna,eleven.int.rna,one.int.rna,three.int.rna)

ys1.1 <- ggplot(odRNAy.int, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
ys1.1

meanrRNAy <- c(mean(odRNAy.int$rRNA[which(odRNAy.int$group_name == 1)]),
                     mean(odRNAy.int$rRNA[which(odRNAy.int$group_name == 2)]),
                     mean(odRNAy.int$rRNA[which(odRNAy.int$group_name == 3)]),
                     mean(odRNAy.int$rRNA[which(odRNAy.int$group_name == 4)]))
GpHy <- GpHy_data$gph
summary(lm(meanrRNAy~GpHy))

ys1.2 <- ggplot(odRNAy.int, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  #geom_jitter(aes(color = ploidy), alpha = 0.5, size = 0.05) +
  geom_jitter(color = "#9b870c", alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  #stat_summary(aes(color = ploidy),fun = mean, geom = "point", size = 4, shape = 21) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  #scale_color_manual(values = c("#9b870c","#643B9F")) +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
ys1.2

summary(lm(odRNAy.int$rRNA~odRNAy.int$GpH))

ys1.3 <- ggplot(odRNAy, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
ys1.3

meanrRNA.ratioy <- c(mean(odRNAy$mRNA.rRNA[which(odRNAy$group_name == "0.92")]),
               mean(odRNAy$mRNA.rRNA[which(odRNAy$group_name == "2.66")]),
               mean(odRNAy$mRNA.rRNA[which(odRNAy$group_name == "5.74")]),
               mean(odRNAy$mRNA.rRNA[which(odRNAy$group_name == "8.61")]))
GpHy <- GpHy_data$gph
summary(lm(meanrRNA.ratioy~GpHy))

ys1.4 <- ggplot(odRNAy, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  geom_jitter(color = "#9b870c", alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
ys1.4

summary(lm(odRNAy$mRNA.rRNA~odRNAy$GpH))

nine.nr.int.rna <- assign_RNA_types_y(nine.nr.int,1,1)
eleven.nr.int.rna <- assign_RNA_types_y(eleven.nr.int,2,2)
one.nr.int.rna <- assign_RNA_types_y(one.nr.int,3,3)
three.nr.int.rna <- assign_RNA_types_y(three.nr.int,4,4)

odRNAy.nr.int <- rbind(nine.nr.int.rna,eleven.nr.int.rna,one.nr.int.rna,three.nr.int.rna)

ys1.5 <- ggplot(odRNAy.nr.int, aes(x = GpH, y = rProt, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
ys1.5

meanrProty <- c(mean(odRNAy.nr.int$rProt[which(odRNAy.nr.int$group_name == 1)]),
                mean(odRNAy.nr.int$rProt[which(odRNAy.nr.int$group_name == 2)]),
                mean(odRNAy.nr.int$rProt[which(odRNAy.nr.int$group_name == 3)]),
                mean(odRNAy.nr.int$rProt[which(odRNAy.nr.int$group_name == 4)]))
GpHy <- GpHy_data$gph
summary(lm(meanrProty~GpHy))

ys1.6 <- ggplot(odRNAy.nr.int, aes(x = GpH, y = rProt, fill = factor(GpH))) +
  geom_jitter(color = "#9b870c", alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
ys1.6

summary(lm(odRNAy.nr.int$rProt~odRNAy.nr.int$GpH))

#bacteria

growcurvb.int <- scale_normalize_seurat_data(growcurvb, 3000, 20, NULL)

OD0.5.int <- scale_normalize_seurat_data(OD0.5, 3000, 20, NULL)
OD1.0.int <- scale_normalize_seurat_data(OD1.0, 3000, 20, NULL)
OD1.3.int <- scale_normalize_seurat_data(OD1.3, 3000, 20, NULL)
OD1.6.int <- scale_normalize_seurat_data(OD1.6, 3000, 20, NULL)
OD2.8.int <- scale_normalize_seurat_data(OD2.8, 3000, 20, NULL)
OD3.6.int <- scale_normalize_seurat_data(OD3.6, 3000, 20, NULL)
OD5.3.int <- scale_normalize_seurat_data(OD5.3, 3000, 20, NULL)
OD6.0.int <- scale_normalize_seurat_data(OD6.0, 3000, 20, NULL)

OD0.5.RNA.int <- assign_RNA_types_b(OD0.5.int, "0.5", 4)
OD1.0.RNA.int <- assign_RNA_types_b(OD1.0.int, "1.0", 5)
OD1.3.RNA.int <- assign_RNA_types_b(OD1.3.int, "1.3", 8)
OD1.6.RNA.int <- assign_RNA_types_b(OD1.6.int, "1.6", 9)
OD2.8.RNA.int <- assign_RNA_types_b(OD2.8.int, "2.8", 12)
OD3.6.RNA.int <- assign_RNA_types_b(OD3.6.int, "3.6", 13)
OD5.3.RNA.int <- assign_RNA_types_b(OD5.3.int, "5.3", 14)
OD6.0.RNA.int <- assign_RNA_types_b(OD6.0.int, "6.0", 15)

odRNAb.int <- rbind(OD0.5.RNA.int,OD1.0.RNA.int,OD1.3.RNA.int,OD1.6.RNA.int,OD2.8.RNA.int,OD3.6.RNA.int,OD5.3.RNA.int,OD6.0.RNA.int)

bs1.1 <- ggplot(odRNAb.int, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.1

x <- c(mean(OD0.5.RNA.int$GpH), 
       mean(OD1.0.RNA.int$GpH), 
       mean(OD1.3.RNA.int$GpH), 
       mean(OD1.6.RNA.int$GpH),
       mean(OD2.8.RNA.int$GpH),
       mean(OD3.6.RNA.int$GpH),
       mean(OD5.3.RNA.int$GpH),
       mean(OD6.0.RNA.int$GpH))

y <- c(mean(OD0.5.RNA.int$rRNA), 
       mean(OD1.0.RNA.int$rRNA), 
       mean(OD1.3.RNA.int$rRNA), 
       mean(OD1.6.RNA.int$rRNA),
       mean(OD2.8.RNA.int$rRNA),
       mean(OD3.6.RNA.int$rRNA),
       mean(OD5.3.RNA.int$rRNA),
       mean(OD6.0.RNA.int$rRNA))

summary(lm(y~x))

bs1.2 <- ggplot(odRNAb.int, aes(x = GpH, y = rRNA, fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.2

summary(lm(odRNAb.int$rRNA~odRNAb.int$GpH))

bs1.3 <- ggplot(odRNAb, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.3

x <- c(mean(OD0.5.RNA$GpH), 
       mean(OD1.0.RNA$GpH), 
       mean(OD1.3.RNA$GpH), 
       mean(OD1.6.RNA$GpH),
       mean(OD2.8.RNA$GpH),
       mean(OD3.6.RNA$GpH),
       mean(OD5.3.RNA$GpH),
       mean(OD6.0.RNA$GpH))

y <- c(mean(OD0.5.RNA$mRNA.rRNA), 
       mean(OD1.0.RNA$mRNA.rRNA), 
       mean(OD1.3.RNA$mRNA.rRNA), 
       mean(OD1.6.RNA$mRNA.rRNA),
       mean(OD2.8.RNA$mRNA.rRNA),
       mean(OD3.6.RNA$mRNA.rRNA),
       mean(OD5.3.RNA$mRNA.rRNA),
       mean(OD6.0.RNA$mRNA.rRNA))

summary(lm(y~x))

bs1.4 <- ggplot(odRNAb, aes(x = GpH, y = mRNA.rRNA, fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.4

summary(lm(odRNAb$mRNA.rRNA~odRNAb$GpH))

OD0.5.nr.int <- scale_normalize_seurat_data(OD0.5.nr, 3000, 20, NULL)
OD1.0.nr.int <- scale_normalize_seurat_data(OD1.0.nr, 3000, 20, NULL)
OD1.3.nr.int <- scale_normalize_seurat_data(OD1.3.nr, 3000, 20, NULL)
OD1.6.nr.int <- scale_normalize_seurat_data(OD1.6.nr, 3000, 20, NULL)
OD2.8.nr.int <- scale_normalize_seurat_data(OD2.8.nr, 3000, 20, NULL)
OD3.6.nr.int <- scale_normalize_seurat_data(OD3.6.nr, 3000, 20, NULL)
OD5.3.nr.int <- scale_normalize_seurat_data(OD5.3.nr, 3000, 20, NULL)
OD6.0.nr.int <- scale_normalize_seurat_data(OD6.0.nr, 3000, 20, NULL)

OD0.5.RNA.nr.int <- assign_RNA_types_b(OD0.5.nr.int, "0.5", 4)
OD1.0.RNA.nr.int <- assign_RNA_types_b(OD1.0.nr.int, "1.0", 5)
OD1.3.RNA.nr.int <- assign_RNA_types_b(OD1.3.nr.int, "1.3", 8)
OD1.6.RNA.nr.int <- assign_RNA_types_b(OD1.6.nr.int, "1.6", 9)
OD2.8.RNA.nr.int <- assign_RNA_types_b(OD2.8.nr.int, "2.8", 12)
OD3.6.RNA.nr.int <- assign_RNA_types_b(OD3.6.nr.int, "3.6", 13)
OD5.3.RNA.nr.int <- assign_RNA_types_b(OD5.3.nr.int, "5.3", 14)
OD6.0.RNA.nr.int <- assign_RNA_types_b(OD6.0.nr.int, "6.0", 15)

odRNAb.nr.int <- rbind(OD0.5.RNA.nr.int,
                       OD1.0.RNA.nr.int,
                       OD1.3.RNA.nr.int,
                       OD1.6.RNA.nr.int,
                       OD2.8.RNA.nr.int,
                       OD3.6.RNA.nr.int,
                       OD5.3.RNA.nr.int,
                       OD6.0.RNA.nr.int)

bs1.5 <- ggplot(odRNAb.nr.int, aes(x = GpH, y = rProt, fill = factor(GpH))) +
  stat_summary(fun = mean, geom = "point", size = 4, fill = "#8CA2CB", shape = 21) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("mean rRNA counts/cell") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.5

x <- c(mean(OD0.5.RNA.nr.int$GpH), 
       mean(OD1.0.RNA.nr.int$GpH), 
       mean(OD1.3.RNA.nr.int$GpH), 
       mean(OD1.6.RNA.nr.int$GpH),
       mean(OD2.8.RNA.nr.int$GpH),
       mean(OD3.6.RNA.nr.int$GpH),
       mean(OD5.3.RNA.nr.int$GpH),
       mean(OD6.0.RNA.nr.int$GpH))

y <- c(mean(OD0.5.RNA.nr.int$rProt), 
       mean(OD1.0.RNA.nr.int$rProt), 
       mean(OD1.3.RNA.nr.int$rProt), 
       mean(OD1.6.RNA.nr.int$rProt),
       mean(OD2.8.RNA.nr.int$rProt),
       mean(OD3.6.RNA.nr.int$rProt),
       mean(OD5.3.RNA.nr.int$rProt),
       mean(OD6.0.RNA.nr.int$rProt))

summary(lm(y~x))

bs1.6 <- ggplot(odRNAb.nr.int, aes(x = GpH, y = rProt, fill = factor(GpH))) +
  geom_jitter(aes(color = "#2c4470"), alpha = 0.5, size = 0.01) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_fill_manual(values = c("#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB","#8CA2CB")) +
  scale_color_manual(values = c("#2c4470")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  NoLegend()
bs1.6

summary(lm(odRNAb.nr.int$rProt~odRNAb.nr.int$GpH))

#--SFs yeast growth/stress, cell cycle, ploidy----------------------------------

growcurvy.nr.int <- scale_normalize_seurat_data(growcurvy.nr, 3000, 20, NULL)
growcurvy.nr.clust <- find_neighbors_clusters(growcurvy.nr.int, 10, 0.2)
nine.nr.int <- scale_normalize_seurat_data(nine.nr, 3000, 20, NULL)
nine.nr.clust <- find_neighbors_clusters(nine.nr.int, 10, 0.2)
eleven.nr.int <- scale_normalize_seurat_data(eleven.nr, 3000, 20, NULL)
eleven.nr.clust <- find_neighbors_clusters(eleven.nr.int, 10, 0.2)
one.nr.int <- scale_normalize_seurat_data(one.nr, 3000, 20, NULL)
one.nr.clust <- find_neighbors_clusters(one.nr.int, 10, 0.2)
three.nr.int <- scale_normalize_seurat_data(three.nr, 3000, 20, NULL)
three.nr.clust <- find_neighbors_clusters(three.nr.int, 10, 0.1)

gasch_spellman_class <- function(SO, SO1) {
  
  ribosomes <- read.csv('C:/Users/lbrettne/Desktop/Solo.out/ribosomes.csv')
  names(ribosomes)[1] = "name1"
  rRNA.genes <- ribosomes[which(ribosomes$type == 'rRNA'),]
  rRNA  <- colSums(SO1@assays$RNA[which(rownames(SO1@assays$RNA) %in% rRNA.genes$name1 | rownames(SO1@assays$RNA) %in% rRNA.genes$name2),])
  totalRNA <- colSums(SO1@assays$RNA)
  
  class <- read.csv('C:/Users/lbrettne/Desktop/RiBi_iESR_CC_genes.csv')
  class$Gene <- convert_gene_names_yeast(class$Gene, FALSE)
  RP.genes <- class[which(class$Group == 'RP'),]
  iESR.genes <- class[which(class$Group == 'iESR'),]
  RiBi.genes <- class[which(class$Group == 'RiBi'),]
  M_G1.genes <- class[which(class$Group == 'M/G1'),]
  G1.genes <- class[which(class$Group == 'G1'),]
  S.genes <- class[which(class$Group == 'S'),]
  G2.genes <- class[which(class$Group == 'G2'),]
  M.genes <- class[which(class$Group == 'M'),]
  Mito.genes <- class[which(class$Group == 'Mito'),]
  
  
  RP  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% RP.genes$Gene),])
  iESR  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% iESR.genes$Gene),])
  RiBi  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% RiBi.genes$Gene),])
  M_G1  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% M_G1.genes$Gene),])
  G1  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% G1.genes$Gene),])
  S  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% S.genes$Gene),])
  G2  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% G2.genes$Gene),])
  M  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% M.genes$Gene),])
  totalmRNA <- colSums(SO@assays$RNA)


  temp <- data.frame(rRNA = rRNA,
                     #rRNA = rRNA/totalRNA,
                     RP = RP/totalmRNA,
                     iESR = iESR,
                     iESR.ratio = iESR/totalmRNA,
                     RiBi = RiBi/totalmRNA,
                     M_G1 = M_G1/totalmRNA,
                     G1 = G1/totalmRNA,
                     S = S/totalmRNA,
                     G2 = G2/totalmRNA,
                     M = M/totalmRNA,
                     timepoint = SO$hour)
  return(temp)
}
gsc_SO <- function(SO, gsc){
  SO$rRNA = log10(gsc$rRNA)
  SO$RP = gsc$RP
  SO$RiBi = gsc$RiBi
  SO$iESR = gsc$iESR
  SO$M_G1 = gsc$M_G1
  SO$G1 = gsc$G1
  SO$S = gsc$S
  SO$G2 = gsc$G2
  SO$M = gsc$M
  
  return(SO)
}

growcurvy.gsc = gasch_spellman_class(growcurvy.nr.int, growcurvy)
growcurvy.nr.clust <- gsc_SO(growcurvy.nr.clust, growcurvy.gsc)
nine.gsc = gasch_spellman_class(nine.nr.clust, nine)
nine.nr.clust <- gsc_SO(nine.nr.clust, nine.gsc)
eleven.gsc = gasch_spellman_class(eleven.nr, eleven)
eleven.nr.clust <- gsc_SO(eleven.nr.clust, eleven.gsc)
one.gsc = gasch_spellman_class(one.nr, one)
one.nr.clust <- gsc_SO(one.nr.clust, one.gsc)
three.gsc = gasch_spellman_class(three.nr, three)
three.nr.clust <- gsc_SO(three.nr.clust, three.gsc)

#Gasch et al. 2017 RP/RiBi/iESR classifications

#full growth curve
DimPlot(growcurvy.nr.clust, reduction = "umap", group.by = "hour")+ 
  theme(aspect.ratio = 1) + scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) +NoLegend()
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "rRNA", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "RP", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "RiBi", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)

#timepoints rRNA
FeaturePlot(nine.nr.clust, reduction = "umap", features = "rRNA", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(eleven.nr.clust, reduction = "umap", features = "rRNA", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(one.nr.clust, reduction = "umap", features = "rRNA", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(three.nr.clust, reduction = "umap", features = "rRNA", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)

#timepoints RP
FeaturePlot(nine.nr.clust, reduction = "umap", features = "RP", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(eleven.nr.clust, reduction = "umap", features = "RP", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(one.nr.clust, reduction = "umap", features = "RP", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(three.nr.clust, reduction = "umap", features = "RP", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)

#timepoints RiBi
FeaturePlot(nine.nr.clust, reduction = "umap", features = "RiBi", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(eleven.nr.clust, reduction = "umap", features = "RiBi", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(one.nr.clust, reduction = "umap", features = "RiBi", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(three.nr.clust, reduction = "umap", features = "RiBi", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(three.nr.clust, reduction = "umap", features = "RiBi")

#timepoints iESR
FeaturePlot(nine.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(eleven.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(one.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(three.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)

#growth curve cell cycle
DimPlot(growcurvy.nr.clust, reduction = "umap", group.by = "hour")+ 
  theme(aspect.ratio = 1) + scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) +NoLegend()
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "rRNA", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "M_G1", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "G1", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "S", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "G2", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "M", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)

summary(lm(log10(growcurvy.gsc$rRNA)~growcurvy.gsc$M_G1+growcurvy.gsc$G1+growcurvy.gsc$S+growcurvy.gsc$G2+growcurvy.gsc$M))
summary(lm(log10(growcurvy.gsc$rRNA)~growcurvy.gsc$M_G1))
summary(lm(log10(growcurvy.gsc$rRNA)~growcurvy.gsc$G1))
summary(lm(log10(growcurvy.gsc$rRNA)~growcurvy.gsc$S))
summary(lm(log10(growcurvy.gsc$rRNA)~growcurvy.gsc$G2))
summary(lm(log10(growcurvy.gsc$rRNA)~growcurvy.gsc$M))

ggplot(growcurvy.gsc, aes(x = G1, y = log10(rRNA))) +
  geom_point(aes(color = timepoint), size = 0.25, alpha = 0.5) +
  scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.75) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()

ggplot(growcurvy.gsc, aes(x = S, y = log10(rRNA))) +
  geom_point(aes(color = timepoint), size = 0.25, alpha = 0.5) +
  scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.75) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()

ggplot(growcurvy.gsc, aes(x = G2, y = log10(rRNA))) +
  geom_point(aes(color = timepoint), size = 0.25, alpha = 0.5) +
  scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.75) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()

ggplot(growcurvy.gsc, aes(x = M, y = log10(rRNA))) +
  geom_point(aes(color = timepoint), size = 0.25, alpha = 0.5) +
  scale_color_manual(values = c("#8A3324","#ED820E","#006666","#90E0EF")) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.75) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()

#growth curve ploidy
growcurv.nr.clust$ploidy = 0
growcurv.nr.clust$ploidy[which(growcurv.nr.clust$celltype == "s288c")] = 1
summary(lm(growcurv.gsc$rRNA~growcurv.nr.clust$ploidy))

DimPlot(growcurvy.nr.clust, reduction = "umap", group.by = "celltype")+ 
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) + scale_color_manual(values = c("#9b870c","#643B9F")) + NoLegend() 

diploids <- odRNAy[which(odRNAy$ploidy == "s288c"),]
haploids <- odRNAy[which(odRNAy$ploidy == "by4741"),]

ggplot(odRNAy, aes(x = GpH, y = log10(rRNA), fill = factor(GpH))) +
  geom_jitter(aes(color = ploidy), alpha = 0.5, size = 0.05) +
  geom_violin(trim = FALSE, alpha = 0) +
  stat_summary(data = diploids, aes(x = GpH, y = log10(rRNA), fill = ploidy),fun = mean, geom = "point", size = 4, shape = 21) +
  stat_summary(data = haploids, aes(x = GpH, y = log10(rRNA), fill = ploidy),fun = mean, geom = "point", size = 4, shape = 21) +
  geom_smooth(data = diploids, aes(x = GpH, y = log10(rRNA), fill = ploidy, color = ploidy), method = "lm", se = FALSE, size = 1.5) +
  geom_smooth(data = haploids, aes(x = GpH, y = log10(rRNA), fill = ploidy, color = ploidy), method = "lm", se = FALSE, size = 1.5, linetype = "dashed") +
  ylab("log10(rRNA counts/cell)") +
  xlab("growth rate (generations/hour)") +
  scale_color_manual(values = c("#9b870c","#643B9F")) +
  scale_fill_manual(values = c("#FFEBC0","#FFEBC0","#FFEBC0","#FFEBC0","#9b870c","#643B9F")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()




#--SFs yeast additional Brauer analyses-----------------------------------------

ys2.1 <- ggplot(all.br, aes(x = GpH, y = bulk, fill = factor(GpH))) +
  #geom_jitter(color = "#99CCB2", alpha = 0.5, size = 0.01) +
  #geom_violin(trim = FALSE, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(aes(fill = GpH), method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("inferred scGR") +
  xlab("PGR (gen/hour)") +
  scale_fill_manual(values = c("#89BB87","#89BB87","#89BB87","#89BB87")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)
  ) +
  NoLegend()
ys2.1

summary(lm(all.br$gr~all.br$GpH))

ys2.2 <- ggplot(baseline_nine, aes(x = gr_norm, y = log10(rRNA))) +
  geom_point(color = "#89BB87") +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 1) +
  ylab("rRNA/cell") +
  xlab("inferred scGR") +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()
  ) +
  scale_y_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  NoLegend()
ys2.2

summary(lm(all.br$rRNA.ratio~all.br$gr))

ys2.1 + ys2.2

retrotransposons <- function(SO, SO1) {
  
  ribosomes <- read.csv('C:/Users/lbrettne/Desktop/Solo.out/ribosomes.csv')
  names(ribosomes)[1] = "name1"
  rRNA.genes <- ribosomes[which(ribosomes$type == 'rRNA'),]
  rRNA  <- colSums(SO1@assays$RNA[which(rownames(SO1@assays$RNA) %in% rRNA.genes$name1 | rownames(SO1@assays$RNA) %in% rRNA.genes$name2),])
  totalRNA <- colSums(SO1@assays$RNA)
  
  rts <- read.csv('C:/Users/lbrettne/Desktop/retrotransposons.csv')
  rts$gene <- convert_gene_names_yeast(rts$gene, FALSE)
  teg.genes <- rts[which(rts$type == 'transposable element gene'),]
  ltr.genes <- rts[which(rts$type == 'LTR retrotransposon'),]
  pseudo.genes <- rts[which(rts$type == 'rt pseudogene'),]
  
  rts  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% rts$gene),])
  teg  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% teg.genes$gene),])
  ltr  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% ltr.genes$gene),])
  pse  <- colSums(SO@assays$RNA[which(rownames(SO@assays$RNA) %in% pseudo.genes$gene),])
  totalmRNA <- colSums(SO@assays$RNA)
  
  temp <- data.frame(rRNA = rRNA,
                     rRNA.ratio = rRNA/totalRNA,
                     rts = rts,
                     rts.ratio = rts/totalmRNA,
                     teg = teg,
                     teg.ratio = teg/totalmRNA,
                     ltr = ltr,
                     ltr.ratio = ltr/totalmRNA,
                     pse = pse,
                     pse.ratio = pse/totalmRNA,
                     timepoint = SO$hour)
  return(temp)
}
rts_SO <- function(SO, rts){
  SO$rRNA = log10(rts$rRNA)
  SO$rRNA.ratio = rts$rRNA.ratio
  SO$rts = rts$rts
  SO$rts.ratio = rts$rts.ratio
  SO$teg = rts$teg
  SO$teg.ratio = rts$teg.ratio
  SO$ltr = rts$ltr
  SO$ltr.ratio = rts$ltr.ratio
  SO$pse = rts$pse
  SO$pse.ratio = rts$pse.ratio
  
  
  return(SO)
}

growcurvy.rts <- retrotransposons(growcurvy.nr.int, growcurvy)
growcurvy.nr.clust <- rts_SO(growcurvy.nr.clust, growcurvy.rts)
nine.rts <- retrotransposons(nine.nr.int, nine)
nine.nr.clust <- rts_SO(nine.nr.clust, nine.rts)
eleven.rts <- retrotransposons(eleven.nr.int, eleven)
eleven.nr.clust <- rts_SO(eleven.nr.clust, eleven.rts)
one.rts <- retrotransposons(one.nr.int, one)
one.nr.clust <- rts_SO(one.nr.clust, one.rts)
three.rts <- retrotransposons(three.nr.int, three)
three.nr.clust <- rts_SO(three.nr.clust, three.rts)

FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(growcurvy.nr.clust, reduction = "umap", features = "rts", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)

FeaturePlot(nine.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(nine.nr.clust, reduction = "umap", features = "rts", cols = c("black","#89BB87", "#FFEBC0"))+ 
  theme(aspect.ratio = 1)

FeaturePlot(eleven.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(eleven.nr.clust, reduction = "umap", features = "rts", cols = c("black","#89BB87", "#FFEBC0"))+ 
  theme(aspect.ratio = 1)

FeaturePlot(one.nr.clust, reduction = "umap", features = "rRNA.ratio", cols = c("black", "purple", "red", "yellow"))+ 
  theme(aspect.ratio = 1)
FeaturePlot(one.nr.clust, reduction = "umap", features = "rts", cols = c("black","#89BB87", "#FFEBC0"))+ 
  theme(aspect.ratio = 1)

FeaturePlot(three.nr.clust, reduction = "umap", features = "iESR", cols = c("black", "yellow"), "red", "yellow"))+ 
  theme(aspect.ratio = 1)

t1 <- FeaturePlot(nine.nr.clust, reduction = "umap", features = "rts", cols = c("black","#89BB87", "#FFEBC0"), pt.size = 1.5)+ 
  theme(aspect.ratio = 1) + NoLegend()
t2 <- FeaturePlot(three.nr.clust, reduction = "umap", features = "rts", cols = c("black","#89BB87", "#FFEBC0"), pt.size = 1.5)+ 
  theme(aspect.ratio = 1)+ NoLegend()
plot_grid(t1,t2,b5.3,b5.4)

DimPlot(three.nr.clust, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1) + 
  NoLegend() + 
  theme(aspect.ratio = 1)

differentially_expressed_genes(three.nr.clust)


length(which(nine.br$gr==0))/length(nine.br$gr)
length(which(eleven.br$gr==0))/length(eleven.br$gr)
length(which(one.br$gr==0))/length(one.br$gr)
length(which(three.br$gr==0))/length(three.br$gr)

percentile <- length(baseline_nine$group_name)*0.25

test <- baseline_nine[order(baseline_nine$gr_norm, decreasing = TRUE),]
test_top <- test[1:percentile,]
test_bottom <- test[(length(test$group_name)-percentile+1):length(test$group_name),]

test_top$percentile <- "top"
test_bottom$percentile <- "bottom"

brauer_25 <- rbind(test_top, test_bottom)

color_top = "#89BB87"
color_bot = "#808080"

ggplot(brauer_25, aes(x = log10(rRNA), fill = percentile)) +
  #geom_histogram(position = "dodge", binwidth = 0.03) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean(log10(test_top$rRNA))), color = color_top, linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = mean(log10(test_bottom$rRNA))), color = color_bot, linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)) +
  scale_fill_manual(values = c(color_bot,color_top)) +
  xlab("rRNA counts/cell") +
  NoLegend()


percentile <- length(baseline_nine$group_name)*0.25
test <- baseline_nine[order(baseline_nine$rsect, decreasing = TRUE),]
test_top <- test[1:percentile,]
test_bottom <- test[(length(test$group_name)-percentile+1):length(test$group_name),]

test_top$percentile <- "top"
test_bottom$percentile <- "bottom"

brauer_25 <- rbind(test_top, test_bottom)

color_top = "#89BB87"
color_bot = "#808080"

ggplot(brauer_25, aes(x = log10(rRNA), fill = percentile)) +
  #geom_histogram(position = "dodge", binwidth = 0.03) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean(log10(test_top$rRNA))), color = color_top,linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = mean(log10(test_bottom$rRNA))), color = color_bot, linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16)) +
  scale_fill_manual(values = c(color_bot,color_top)) +
  xlab("rRNA counts/cell") +
  NoLegend()
  

percentile <- ceiling(length(odRNAb$group_name)*0.25)
test <- odRNAb[order(odRNAb$rsect, decreasing = TRUE),]
test_top <- test[1:percentile,]
test_bottom <- test[(length(test$group_name)-percentile+1):length(test$group_name),]

test_top$percentile <- "top"
test_bottom$percentile <- "bottom"

brauer_25 <- rbind(test_top, test_bottom)

color_top = "#89BB87"
color_bot = "#808080"

ggplot(brauer_25, aes(x = log10(rRNA), fill = percentile)) +
  #geom_histogram(position = "dodge", binwidth = 0.03) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = mean(log10(test_top$rRNA))), color = color_top,linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = mean(log10(test_bottom$rRNA))), color = color_bot, linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),) +
  scale_fill_manual(values = c(color_bot,color_top)) +
  xlab("rRNA counts/cell") +
  NoLegend()

rep1$GpH <- as.factor(rep1$GpH)
rep2$GpH <- as.factor(rep2$GpH)
levels(rep1$GpH) <- c("0.09", "0.35", "0.78", "1.05")  
levels(rep2$GpH) <- c("0.09", "0.35", "0.78", "1.05")  


ggplot(rep2, aes(x = log10(rRNA), y = as.factor(GpH))) + 
  stat_density_ridges(fill = "#9b870c", alpha = 0.5, color = "black", size = 1) + 
  stat_density_ridges(data = rep1, aes(x = log10(rRNA), y = as.factor(GpH)), fill = "#FFEBC0", alpha = 0.5, color = "#9b870c", size = 1) +
  scale_x_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()) +
  scale_fill_manual(values = c(color_bot,color_top)) +
  xlab("rRNA counts/cell") +
  ylab("growth rate (generations/hour)")

odRNAb$GpH <- as.factor(odRNAb$GpH)
levels(odRNAb$GpH) <- c("0.16", "0.33", "0.69", "0.85", "1.07", "1.14", "1.24", "1.27")
odRNAb.m14$GpH <- as.factor(odRNAb.m14$GpH)
levels(odRNAb.m14$GpH) <- c("0.76", "0.85", "0.97", "1.00", "1.24", "1.27")

ggplot(odRNAb.m14, aes(x = log10(rRNA), y = as.factor(GpH))) + 
  stat_density_ridges(fill = "#2c4470", color = "black", alpha = 0.5, size = 1) + 
  stat_density_ridges(data = odRNAb, aes(x = log10(rRNA), y = as.factor(GpH)), fill = "#8CA2CB", color = "#2c4470", alpha = 0.5, size = 1) +
  scale_x_continuous(breaks = c(3,4,5),labels = c("1,000", "10,000", "100,000")) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16),
    axis.title = element_blank()) +
  scale_fill_manual(values = c(color_bot,color_top)) +
  xlab("rRNA counts/cell") +
  ylab("growth rate (generations/hour)")


test <- data.frame(x = rbeta(100000, 1.2, 5))

ggplot(test, aes(x = x)) +
  geom_density(alpha = 0.5, fill = "#89BB87") +
  geom_vline(aes(xintercept = mean(x)), color = "black",linetype = "dashed", size = 1) +
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    axis.line = element_line(size = 0.5, color = "black"),
    aspect.ratio = 1,
    text = element_text(size=16))
