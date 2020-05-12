---
#title: "R Script for WGCNA"
#wd: choose your own working directory, and put the path into setwd(" ") at 12th line
#input: all.gene_counts.tsv (raw counts of RNAseq from 36 samples)
#output: 
---

#########################################################
# WGCNA analysis 
# moderate cut-off, signed model, dicor, protein_coding gene

setwd("./") # set working directory

# attach required packages
library(WGCNA)
library(DESeq2)
library(sva)
library(factoextra)
library(reshape2)
library(circlize)
library(limma)
library(org.Hs.eg.db)
library(edgeR)
library(ComplexHeatmap)
library(gplots)
options(stringsAsFactors = FALSE)

#=================================================
# 1. Data input and cleaning 
## read in raw counts file 
all_counts <- read.delim2("./all.gene_counts.tsv",stringsAsFactors=FALSE)
all_gene_counts <- all_counts[(!duplicated(all_counts$ensembl_gene_id)),] # remove columns with duplicated ensembl_gene_id
all_gene_counts <- all_gene_counts[,c(1,8:43)] # remove the annotation columns, only keep gene expression matrix
rownames(all_gene_counts) <- all_gene_counts[,1]
all_gene_counts <- all_gene_counts[,-1]
colnames(all_gene_counts) <- gsub('sample.','',colnames(all_gene_counts))
gene_info <- all_counts[,c(1:4,7)]
rm(all_counts)

##  Normalize: log-transform cpm (cpm , log2(x+1))
all_gene_cpm <- cpm(all_gene_counts)
all_gene_cpm <- log2(all_gene_cpm + 1)

# WGCNA step1 function 
WGCNA_step1 <- function(wd, log2_cpm_df, counts_df, gene_info_df, low_count = 0, protein_coding){
  ## Input
  # wd: working directory
  # log2_cpm_df: gene expression df, log2(cpm + 1), row: gene, column: samples
  # count_df: gene expression df, raw counts
  # gene_info_df: gene information df, contains ensembl_gene_id, gene_biotype(protein coding ...)
  # low_count: int, threshold to filter low counts genes
  # protein_coding: 'protein_coding', 'non_coding', else (keep protein_coding, keep non_coding, no filter)
  
  ## Output
  # two dendrogram figure (before/after remove 1 sample)
  # save 2 version of RData, log2_cpm_df, gene_info_df, sample_info_repeat_order
  
  ## Filter low counts
  # counts <10 in >90% samples
  keep_gene <- rowSums(counts_df < low_count) <= 0.9*ncol(counts_df) 
  print('filter of low counts:')
  print(table(keep_gene))
  log2_cpm_df <- log2_cpm_df[keep_gene,]
  
  ## keep protwin genes / non_coding genes
  if (protein_coding == 'protein_coding'){
    # protein codings 
    keep_protein_gene <- rownames(log2_cpm_df) %in% gene_info_df$ensembl_gene_id[gene_info_df$gene_biotype == 'protein_coding']
    print('romove non-coding genes:')
    print(table(keep_protein_gene))
    log2_cpm_df <- log2_cpm_df[keep_protein_gene,] 
  } else if (protein_coding == 'non_coding'){
    # non coding genes 
    keep_noncoding_gene <- rownames(log2_cpm_df) %in% gene_info_df$ensembl_gene_id[gene_info_df$gene_biotype != 'protein_coding']
    print('remove protein-coding genes:')
    print(table(keep_noncoding_gene))
    log2_cpm_df <- log2_cpm_df[keep_noncoding_gene,] 
  }
  
  ## Transpose expression data
  log2_cpm_df0 <- as.data.frame(t(log2_cpm_df))
  names(log2_cpm_df0) = rownames(log2_cpm_df)
  rownames(log2_cpm_df0) = colnames(log2_cpm_df)
  
  ## check for missing data
  gsg = goodSamplesGenes(log2_cpm_df0, verbose = 3) 
  print('gsg$allOK')
  print(gsg$allOK)
  
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    log2_cpm_df0 = log2_cpm_df0[gsg$goodSamples, gsg$goodGenes]
  }
  
  ## cluster the samples & check for outliers
  sampleTree = hclust(dist(log2_cpm_df0), method = "average")
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  pdf('./plot/sample_cluster_info_all_dendrogram.pdf',width = 7.25, height = 4.64)
  #sizeGrWindow(12,9)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  dev.off()
  
  ## Remove 33947_1
  log2_cpm_df <- log2_cpm_df0[-match('33947_1',rownames(log2_cpm_df0)),] # 35 samples
  
  ##  Loading clinical trait data
  sample_info <- data.frame("sample" = c('4831','4861','4855','4717','4829','4857',
                                         '30013','4230','2147','4198','2173','33947'),
                            "condition" = c(rep('pre_onset',6),rep('post_onset',6)),
                            "gender" = c(rep('male',3),rep('female',3),rep('male',3),rep('female',3)),
                            "age" = c(23,16,11,44,21,23,54,55,55,63,52,71))
  sample_info_repeat <- sample_info[rep(row.names(sample_info),each = 3),]
  rownames(sample_info_repeat) <- paste0(sample_info_repeat$sample,'_',c(1,2,3))
  # change order of sample_info_repeat
  sample_info_repeat_order <- sample_info_repeat[rownames(log2_cpm_df),]
  
  ## show relation between sample and traits
  # Re-cluster samples
  sampleTree = hclust(dist(log2_cpm_df), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(apply(sample_info_repeat_order[,c(1,4)],2,as.numeric), signed = FALSE)
  traitColors = cbind(traitColors,ifelse(sample_info_repeat_order$condition == 'pre_onset','#4B0055','#FDE333'))
  traitColors = cbind(traitColors,ifelse(sample_info_repeat_order$gender == 'female','#8700F9','#00C4AA'))
  age_color <- colorRamp2(c(0,80), c('white', 'red'))
  traitColors = cbind(traitColors,age_color(sample_info_repeat_order$age))
  traitColors <- traitColors[,c(-1,-2)]
  # Plot the sample dendrogram and the colors underneath.
  pdf('./plot/sample_cluster_info_dendrogram.pdf',width = 7.25, height = 4.64)
  plotDendroAndColors(sampleTree, traitColors,groupLabels = c('condition','gender','age'))
  dev.off()
  
  save(log2_cpm_df, gene_info_df, sample_info_repeat_order, file = "01-dataInput.RData")
  # version=2
  save(log2_cpm_df, gene_info_df, sample_info_repeat_order, file = "01-dataInput_V2.RData",version = 2)
  
}

# non-coding gene
WGCNA_step1("E:/WUSTL/Rotation/Yoo Lab/HD/signed_network_correct/",
            log2_cpm_df = all_gene_cpm,
            counts_df = all_gene_counts,
            gene_info_df = gene_info, 
            low_count = 10, 
            protein_coding = 'protein_coding')


#=================================================
# 2. One step network construction and module detection
# rm(list = ls())
load("E:/WUSTL/Rotation/Yoo Lab/HD/non_coding/01-dataInput.RData") # load the output file of 1st step
# log2_cpm_df,gene_info_df,sample_info_repeat_order 
choose_soft_threshold <- function(log2_cpm_df, power_to, type , corType) {
  ## Input 
  # log2_cpm_df: gene expression df, gene expression df, log2(cpm + 1), row: gene, column: samples
  # power_to: int, Choose a set of soft-thresholding powers, in case there is no parallel backend
  # type = 'signed', 'unsigned', 'signed hybrid'
  # corType = "pearson", "bicor"
  
  ## Output
  # sft
  
  ## Choosing the soft-thresholding power:
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to= power_to, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(log2_cpm_df, powerVector = powers, verbose = 5, networkType = type)
  
  # Plot the results:
  pdf('./plot/pick_soft-threshold.pdf',width = 8.58, height = 5.00)
  #sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  # Input of net using HPC
  nGenes = ncol(log2_cpm_df)
  nSamples = nrow(log2_cpm_df)
  maxPOutliers = ifelse(corType=="pearson",1,0.05) # 0.05
  
  save(sft, nGenes, nSamples, type, corType, maxPOutliers, file = 'net_input.RData', version = 2)
}

choose_soft_threshold(log2_cpm_df = log2_cpm_df, power_to = 20,
                      type = 'signed', corType = "bicor")

## One-step network construction and module detection (HPC)
# power = sft$powerEstimate 
### calculate using HPC ->
#### /casa/labs/lab_yoo/users/shuangjia/Project/HD/non_coding_signed
# HD_signed_TOM-block.1.RData
# 
#### /casa/labs/lab_yoo/users/shuangjia/Project/HD/signed_network/net.RData
net = blockwiseModules(log2_cpm_df, power = sft$powerEstimate,
                       networkType = type,
                       TOMType = 'unsigned', minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, corType = corType,
                       saveTOMFileBase = "HD_signed_TOM",
                       verbose = 3,
                       maxBlockSize = nGenes,
                       maxPOutliers = maxPOutliers) 
save(net, file = 'net.RData')


load("./HPC_result/net.RData") # load RData from HPC
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 9.39 * 4.44 # ./plot/gene_clustering_dendrogram.pdf


# save the module assignment and module eigengene information
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "02-networkConstruction-auto.RData")

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "02-networkConstruction-auto_V2.RData", version = 2)
# version2 is used to run in this server: lab_yoo@regmedhpc2.wustl.edu

#=================================================
# 3. Relating modules to external clinical traits
load("01-dataInput.RData")
load("02-networkConstruction-auto.RData")
load("./HPC_result/net.RData")
load("net_input.RData")

trait_char_to_int <- function(sample_info_repeat_order){
  # transfer trait matrix from character to int
  sample_info_int <- sample_info_repeat_order
  sample_info_int$sample <- as.numeric(sample_info_int$sample)
  # pre_oneset - 0; post_onset - 1
  sample_info_int$condition <- ifelse(sample_info_int$condition == 'pre_onset',0,1)
  # female - 0; male - 1
  sample_info_int$gender <- ifelse(sample_info_int$gender == 'female', 0, 1)
  return(sample_info_int)
}
sample_info_int <- trait_char_to_int(sample_info_repeat_order)

correlation_modules_traits <- function(log2_cpm_df,sample_info_int){
  # Define numbers of genes and samples
  nGenes = ncol(log2_cpm_df)
  nSamples = nrow(log2_cpm_df)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(log2_cpm_df, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  
  ## Calculate correlation between modules and traits
  # moduleTraitCor = cor(MEs, sample_info_int, use = "p")
  # moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  corType = "bicor"
  robustY = ifelse(corType=="pearson",T,F)
  if (corType=="pearson") {
    modTraitCor = cor(MEs, traitData, use = "p")
    modTraitP = corPvalueStudent(modTraitCor, nSamples)
  } else {
    modTraitCorP = bicorAndPvalue(MEs, sample_info_int, robustY=robustY)
    modTraitCor = modTraitCorP$bicor
    modTraitP   = modTraitCorP$p
  }
  
  # Plot of correlation between modules and traits
  pdf('./plot/module_trait_relationships.pdf',width = 14.28, height = 3.9)
  # Will display correlations and their p-values
  textMatrix = paste(signif(modTraitCor, 2), "\n(",
                     signif(modTraitP, 1), ")", sep = "");
  dim(textMatrix) = dim(modTraitCor)
  par(mar = c(6, 8.5, 3, 3))
  moduleTraitCor_t <- t(modTraitCor)
  textMatrix_t <- t(textMatrix)
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor_t,
                 yLabels = names(sample_info_int),
                 xLabels = names(MEs),
                 xSymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix_t,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  # 14.28* 3.9 #./plot/WGCNA/module_trait_relationships.pdf
  dev.off()
  save(modTraitCor, modTraitP, file = '03-modules.trait.genes-correlation_00.RData')
}
correlation_modules_traits(log2_cpm_df,sample_info_int)

correlation_gene_traits <- function(){
  ## Calculate Gene relationship to trait and important modules
  # Define variable weight containing the weight column of datTrait
  condition = as.data.frame(sample_info_int$condition)
  names(condition) = "condition"
  
  MEs0 = moduleEigengenes(log2_cpm_df, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  
  robustY = ifelse(corType=="pearson",T,F)
  
  # Gene & Modules
  if (corType=="pearsoon") {
    geneModuleMembership = as.data.frame(cor(log2_cpm_df, MEs, use = "p"))
    MMPvalue = as.data.frame(corPvalueStudent(
      as.matrix(geneModuleMembership), nSamples))
  } else {
    geneModuleMembershipA = bicorAndPvalue(log2_cpm_df, MEs, robustY=robustY)
    geneModuleMembership = geneModuleMembershipA$bicor
    MMPvalue   = geneModuleMembershipA$p
  }
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  
  # Gene & Trait
  if (corType=="pearsoon") {
    geneTraitSignificance = as.data.frame(cor(log2_cpm_df, condition, use = "p"))
    GSPvalue = as.data.frame(corPvalueStudent(
      as.matrix(geneTraitCor), nSamples))
  } else {
    geneTraitCorA = bicorAndPvalue(log2_cpm_df, condition, robustY=robustY)
    geneTraitSignificance = as.data.frame(geneTraitCorA$bicor)
    GSPvalue   = as.data.frame(geneTraitCorA$p)
  }
  names(geneTraitSignificance) = paste("GS.", names(condition), sep="")
  names(GSPvalue) = paste("p.GS.", names(condition), sep="")
  
  save(geneModuleMembership, MMPvalue, geneTraitSignificance, GSPvalue,
       MEs, condition,modNames, file = '03-modules.trait.genes-correlation.RData')
}
correlation_gene_traits()

load("03-modules.trait.genes-correlation.RData")
MM_GS_plot <- function(module, color = module){
  # module: Input module name, should be a color
  #Intramodular analysis: identifying genes with high GS and MM
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  pdf('./plot/Module_membership_VS_gene_significance.pdf',width = 5.71, height = 5.52)
  #sizeGrWindow(7, 7)
  par(mfrow = c(1,1))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for condition",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = color)
  # 5.71* 5.52 # ./plot/Module_membership_VS_gene_significance.pdf
  dev.off()
}
MM_GS_plot(module = "brown")

#names(log2_cpm_df)[moduleColors=="yellow"]

#=================================================
#¡¡5. Visualizing the gene network (time assuming, drawing use HPC)

## Gene network (Calculate using HPC)
# Loading objects: TOM
load("./HPC_result/HD_signed_TOM-block.1.RData")
#load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
pdf("Network_heatmap.pdf", width = 9, height =9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

eigengenes_network <- function(MEs,condition ){
  #  Visualizing the network of eigengenes
  ## Input
  # MEs: module eigengenes
  # condition: Isolate condition from the clinical traits
  
  # Add the condition to existing module eigengenes
  MET = orderMEs(cbind(MEs, condition))
  # Plot the relationships among the eigengenes and the trait
  pdf('./plot/network_of_eigengenes_and_condition.pdf',width = 6, height = 8)
  par(cex = 0.9)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2),
                        cex.lab = 0.8, xLabelsAngle = 90)
  # 6* 8 # ./plot/network_of_eigengenes_and_condition.pdf
  dev.off()
}
eigengenes_network(MEs,condition)


#=================================================
# 6. Exporting to Cytoscape (Calculate use HPC) (have not)

# usually should use HPC!!!
export_to_cytoscape <- function(module){
  ## export one module
  
  # import TOM!!! HPC
  load("./HPC_result/HD_signed_TOM-block.1.RData")
  TOM <- as.matrix(TOM) 
  
  # Select module probes
  probes = names(log2_cpm_df)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  # save.image('./WGCNA_v2.RData', version = 2)
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  rm(TOM)
  
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-",module,"_modules.txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-",module,"_modules.txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
}

export_to_cytoscape('yellow')

#=================================================
## a. barplot for condition (trait) and modules
## b. sample number of each modules - barplot
modTraitCor
modTraitP
sample_info_int
MEs

load("03-modules.trait.genes-correlation_00.RData")

condition_modules_barplot <- function(modTraitCor, modTraitP, moduleColors, cut_line = 0.05){
  ## Input, cut_line: float, cut line for p_value = 0.05
  #change matrix to df
  moduleTraitCor_df <- data.frame(modTraitCor)
  moduleTraitCor_df$color <- substr(rownames(moduleTraitCor_df),3,20)
  moduleTraitCor_df$color <- factor(moduleTraitCor_df$color, levels=unique(moduleTraitCor_df$color))
  moduleTraitCor_df$add_color <- ifelse(data.frame(modTraitP)$condition < cut_line, substr(rownames(moduleTraitCor_df),3,20), 'white')
  # moduleTraitCor_df$add_color_2 <- ifelse(modTraitP < 0.0005, substr(rownames(moduleTraitCor_df),3,20), 'white')
  # add_color - P 0.05; # add_color_2 - P 0.0005
  
  ## a. barplot for condition (trait) and modules
  #pdf('./plot/condition_modules_correlation.pdf', width = 11.4, height = 5.55);
  ggplot(moduleTraitCor_df, aes(x = color, y = condition, fill = color, color = 'black')) +
    geom_bar(stat="identity", width=0.8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.3)) +
    ylim(-1,1) +
    geom_hline(yintercept = c(0.33,-0.33),linetype="dashed", color = "red") +
    scale_fill_manual(values = moduleTraitCor_df$add_color) +
    scale_color_manual(values =  'black');
  # 11.4* 5.55 #  ./plot/condition_modules_correlation.pdf
  ggsave(paste0("./plot/condition_modules_correlation_", as.character(cut_line), ".pdf"), width = 11.4, height = 5.55)
  #dev.off()
  
  # ggplot(moduleTraitCor_df, aes(x = color, y = condition, fill = color, color = 'black')) +
  #   geom_bar(stat="identity", width=0.8) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.3)) +
  #   ylim(-1,1) +
  #   geom_hline(yintercept = c(0.55,-0.55),linetype="dashed", color = "red") +
  #   scale_fill_manual(values = moduleTraitCor_df$add_color_2) +
  #   scale_color_manual(values =  'black')
  # # r^2 = 0.55, P = 0.0005
  # # 11.4* 5.55 #  ./plot/condition_modules_correlation_2.pdf
  
  # b. sample number of each modules - barplot
  ## Input - moduleColors
  
  module_number <- data.frame(table(moduleColors),stringsAsFactors = F)
  module_number$moduleColors <- factor(module_number$moduleColors, levels = moduleTraitCor_df$color)
  
  ggplot(module_number, aes(x = moduleColors, y = Freq, fill = moduleColors, color = 'black')) +
    geom_bar(stat="identity", width=0.8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.3)) +
    scale_fill_manual(values = col2hex(moduleTraitCor_df$color)) +
    scale_color_manual(values =  'black')
  # 11.4* 5.55 #  ./plot/module_genenumber.pdf
  ggsave("./plot/module_genenumber.pdf", width = 11.4, height = 5.55)
}
condition_modules_barplot(modTraitCor, modTraitP, moduleColors)

#=================================================
# GO, KEGG enrichment for yellow, and darkred modules (no need)
#¡¡no GO for non-coding genes

find_go_kegg <- function(module_name, type){
  ## Input: module_name = one color, type = go or kegg
  # Select module probes
  in_Module = is.finite(match(moduleColors, module_name))
  probes = names(log2_cpm_df)
  module_probes = probes[in_Module]
  print(paste0(length(module_probes),' genes in ',module_name,' module'))
  
  # change ENSEMBL to ENTREZID
  module_probes_entrezid = gene_info_df$entrezgene[match(module_probes,gene_info_df$ensembl_gene_id)]
  module_probes_entrezid = module_probes_entrezid[!is.na(module_probes_entrezid)]
  print(paste0(length(module_probes_entrezid),' genes in ',module_name,' module after ENTREZID'))
  
  if(type == 'go'){
    # GO enrichment
    go.module = goana(module_probes_entrezid,species = "Hs")
    go.module$P.adjust <- p.adjust(go.module$P.DE, method = "fdr", n = length(go.module$P.DE))
    result = go.module
  }
  else{
    # KEGG enrichment
    kegg.module <- kegga(module_probes_entrezid,species = "Hs")
    kegg.module$P.adjust <- p.adjust(kegg.module$P.DE, method = "fdr", n = length(kegg.module$P.DE))
    result = kegg.module
  }
  
  n = min(sum(result$P.DE <= 0.05), 10)
  result.order <- result[order(result$P.DE),]
  result.order$Term <- factor(result.order$Term, levels=rev(result.order$Term))
  ggplot(result[1:n,],aes(Term,-log(P.DE))) +
    geom_bar(stat="identity",width=0.8) +
    coord_flip() +
    ggtitle(paste0(module_name," module",' ',type)) +
    theme(axis.text.y = element_text(size = 8,color="black"))
  ggsave(paste0("./plot/",module_name,"_",type,"_barplot.pdf"), width = 6.49, height =4.09)
  return(result)
}

go.yellow <- find_go_kegg('yellow','go')
# 6.49* 4.09 # ./plot/Yello_GO_barplot.pdf
# [1] "440 genes in yellow module"
# [1] "1 genes in yellow module after ENTREZID"
# kegg.yellow <- find_go_kegg('yellow','kegg')


# 4 interested GO terms
ggplot(go.yellow[c("GO:0043065","GO:0043066","GO:0048102","GO:0007049"),],aes(Term,-log(P.DE))) +
  geom_bar(stat="identity",width=0.8) +
  coord_flip() +
  ggtitle("Blue module Go (interested)") +
  theme(axis.text.y = element_text(size = 8,color="black"))
# 6.06 * 2.49 # ./plot/WGCNA/Interested_Go_Blue_module.pdf


## get genes from GO terms (not work yet)
# Get the entrez gene identifiers that are mapped to a GO ID
x <- org.Hs.egGO
mapped_genes <- mappedkeys(x)
# Convert to a list
gene_GO <- as.list(x[mapped_genes])
## get GO terms from genes:
# Convert to a list
GO_gene <- as.list(org.Hs.egGO2EG)
# # For org.Hs.egGO2ALLEGS
xx <- as.list(org.Hs.egGO2ALLEGS)
go_id = GOID( GOTERM[ Term(GOTERM) == "chromatin remodeling"])
genes <- get('GO:1900119',org.Hs.egGO2EG)
unlist(mget(genes,org.Hs.egSYMBOL))

#=================================================
# Gene heatmap in select modules 

gene_heatmap_module <- function(module_vector){
  ## Input module_vector - c('module_1','module_2'...)  (module is color)
  
  ## Calculate cpm again
  ## read in raw counts file 
  all_counts <- read.delim2("E:/WUSTL/Rotation/Yoo Lab/HD/all.gene_counts.tsv",stringsAsFactors=FALSE)
  all_gene_counts <- all_counts[(!duplicated(all_counts$ensembl_gene_id)),]
  all_gene_counts <- all_gene_counts[,c(1,8:43)]
  rownames(all_gene_counts) <- all_gene_counts[,1]
  all_gene_counts <- all_gene_counts[,-1]
  colnames(all_gene_counts) <- gsub('sample.','',colnames(all_gene_counts))
  rm(all_counts)
  
  ##  Normalize: cpm
  all_gene_cpm_ori <- cpm(all_gene_counts) # Important, no log2(+1)
  rm(all_gene_counts)
  
  all_gene_cpm_scale <- apply(all_gene_cpm_ori, 1,scale)
  all_gene_cpm_scale <- t(all_gene_cpm_scale)
  rownames(all_gene_cpm_scale) <- rownames(all_gene_cpm_ori)
  colnames(all_gene_cpm_scale) <- colnames(all_gene_cpm_ori)
  rm(all_gene_cpm_ori)
  # remove 33947_1
  all_gene_cpm_scale <- all_gene_cpm_scale[,-21]
  
  # Select genes in input modules
  probes = names(log2_cpm_df)
  
  # Select module probes
  in_Module = is.finite(match(moduleColors, module_vector))
  in_probes = probes[in_Module] # 3139 gene
  in_Module_gene_cpm_scale <- all_gene_cpm_scale[in_probes,]
  
  # heatmap
  pdf('./plot/heatmap_select_modules.pdf', width = 9.36, height = 6.75)
  set.seed(123)
  sample_info_heatmap <- sample_info_int
  # female -0; male - 1
  # pre_oneset - 0; post_onset - 1
  sample_info_heatmap$condition_name <- ifelse(sample_info_heatmap$condition == 0,'pre_oneset','post_onset')
  col_ha_split <- HeatmapAnnotation(condition = sample_info_heatmap$condition_name,
                                    col = list(condition = c('pre_oneset' = '#4B0055',
                                                             'post_onset' = '#FDE333')))
  # col2hex
  row_ha_split <- rowAnnotation(module = moduleColors[match(rownames(in_Module_gene_cpm_scale),names(moduleLabels))],
                                col = list(module = sapply(module_vector,col2hex)))
  print(Heatmap(in_Module_gene_cpm_scale, 
                cluster_rows = T,clustering_method_rows = 'ward.D', show_row_names = F,
                cluster_columns = F,
                col = colorRamp2(c(-3,0,3),c('#0704FF','#EEEEEE','#FF0905')), 
                right_annotation = row_ha_split,
                gap = unit(1,"mm"),
                row_split = moduleColors[match(rownames(in_Module_gene_cpm_scale),names(moduleLabels))],
                top_annotation = col_ha_split,
                column_split = factor(sample_info_heatmap$condition_name, levels = unique(sample_info_heatmap$condition_name)),
                name = 'scaled_cpm'
  )) 
  dev.off()
  # 9.36 * 6.75 # ./plot/heatmap_select_modules.pdf
  
}
gene_heatmap_module(module_vector = c('yellow','brown','green','red','black','greenyellow'))

#=================================================
# Export gene of each modules (format: gene + foldchange)
moduleColors
gene_info_df

# function to export specific module
export_module <- function(module_name){
  # Recalculate all_gene_cpm_ori
  all_counts <- read.delim2("E:/WUSTL/Rotation/Yoo Lab/HD/all.gene_counts.tsv",stringsAsFactors=FALSE)
  all_gene_counts <- all_counts[(!duplicated(all_counts$ensembl_gene_id)),]
  all_gene_counts <- all_gene_counts[,c(1,8:43)]
  rownames(all_gene_counts) <- all_gene_counts[,1]
  all_gene_counts <- all_gene_counts[,-1]
  colnames(all_gene_counts) <- gsub('sample.','',colnames(all_gene_counts))
  rm(all_counts)
  # cpm
  all_gene_cpm_ori <- cpm(all_gene_counts) # Important, no log2(+1)
  
  options(stringsAsFactors = FALSE)
  probes = names(log2_cpm_df)
  in_Module = is.finite(match(moduleColors, module_name))
  module_probes = probes[in_Module]
  # transform
  module_probes_genename = gene_info_df$external_gene_name[match(module_probes,gene_info_df$ensembl_gene_id)]
  #module_probes_genename = module_probes_genename[!is.na(module_probes_genename)]
  module_gene_ex <- all_gene_cpm_ori[module_probes,-21]
  module_gene_ex_pre <- module_gene_ex[,rownames(sample_info_repeat_order[sample_info_repeat_order$condition == 'pre_onset',])]
  module_gene_ex_post <- module_gene_ex[,rownames(sample_info_repeat_order[sample_info_repeat_order$condition == 'post_onset',])]
  module_df <- data.frame('gene_name' = module_probes_genename,
                          'gene_log2fc' = log2(apply(module_gene_ex_post,1,mean)/apply(module_gene_ex_pre,1,mean)))
  write.table(module_df,paste0(module_name,'_module_df.txt'),sep = '\t',row.names = F,quote = F)
  
}

export_module("greenyellow")
c('yellow','brown','green','red','black','greenyellow')


# gene_module_df <- data.frame('external_gene_name' = probes,
#                              'modulecolors' = moduleColors,
#                              'ensembl_gene_id' = gene_info$ensembl_gene_id[match(probes,gene_info$external_gene_name)],
#                              'entrezgene' = gene_info$entrezgene[match(probes,gene_info$external_gene_name)],
#                              'gene_biotype' = gene_info$gene_biotype[match(probes,gene_info$external_gene_name)],
#                              'description' = gene_info$description[match(probes,gene_info$external_gene_name)])
# write.csv(gene_module_df,'gene_module_df')





