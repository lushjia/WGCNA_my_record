# WGCNA
### WGCNA for RNAseq data
This repositories record my opinioin on how to use WGCNA package and some following analysis of the coexpression network. Most information posted here is gathered from this official [WGCNA website]( https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/index.html).
In that website, *Tutorials* part contains all the R script one needed to install and run WGNCA. *FAQ: Problems installing or using the package* part contains some details of how to choose parameters based on the network you want to build. 

#### Some Terms 
- **soft threshold**: a power used to weighted the correlation, calculated by `pickSoftThreshold` function in WGCNA package
- **module**: expressions of genes in one module are highly correlation with each other. In signed network, the genes in one module are positively correlated, while in unsigned network, the genes are either positively or negatively correlated.
- **module eigengene**: PC1 of genes in one module, which could capture main features of the whole module
- **module membership:**: correlation between a gene and the module eigengene, could be used to identify intramodular hub genes
- **adjacency matrix**: weighted correlation matrix among genes
- **TOM(Topological overlap matrix)**: To minimize effects of noise and spurious associations, they transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity

#### Input and Parameters
- **sample size**: They do not recommend attempting WGCNA on a data set consisting of fewer than 15 samples.  If at all possible, one should have at least 20 samples
- **gene expression matrix**: Rows represent the genes. Columns represent the samples. This matrix will be transposed in the WGCNA script. They then recommend a variance-stabilizing transformation. For example, package DESeq2 implements the function varianceStabilizingTransformation which we have found useful, but one could also start with normalized counts (or RPKM/FPKM data) and log-transform them using log2(x+1). For highly expressed features, the differences between full variance stabilization and a simple log transformation are small.
- **trait matrix**: the values of this matrix must be numbers, but not characters
- 











### WGCNA_functins.R
This is an example of using WGCNA to build a signed network for protein coding genes. In this R script, I integrated each step in WGCNA Tutorials into a function, so that we can run each step using one-line command. 

**Notice**: 
1. some commands in this file targets specific on this file: `yooa@storage1.ris.wustl.edu:/Active/Youngmi/htcf.wustl.edu/files/pd67E8el/Yoo_s4611_MGI0029/all.gene_counts.tsv`. 
If you try to run WGCNA on another file, the first step, *1. Data input and cleaning*, needs to be changed. 
2. when running some step, we'd better compute using high performance cluster(HPC) but not our own laptop, because the memory requirment is super large. I marked those steps using (HPC) in my script, and if one try to run it on HPC, he/she should first upload the outputs to the HPC and then run the script. 
