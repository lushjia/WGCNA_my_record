# WGCNA
### WGCNA for RNAseq data
This repositories introduce my opinioin on how to use WGCNA package and some following analysis of the coexpression network. Most information posted here is from this official [WGCNA website]( https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/index.html).
In that website, *Tutorials* part contains all the R script one needed to install and run WGNCA. *FAQ: Problems installing or using the package* part contains some details of how to choose parameters based on the network you want to build. 

### WGCNA_functins.R
This is an example of running WGCNA. In this R script, I integrated each step in WGCNA Tutorials into a function, so that we can run each step using one-line command. 

**Notice**: 
1. some commands in this file targets specific on this file: `yooa@storage1.ris.wustl.edu:/Active/Youngmi/htcf.wustl.edu/files/pd67E8el/Yoo_s4611_MGI0029/all.gene_counts.tsv`. 
If you try to run WGCNA on another file, the first step, *1. Data input and cleaning*, needs to be changed. 
2. when running some step, we'd better compute using high performance cluster(HPC) but not our own laptop, because the memory requirment is super large. I marked those steps using (HPC) in my script, and if one try to run it on HPC, he/she should first upload the outputs to the HPC and then run the script. 
