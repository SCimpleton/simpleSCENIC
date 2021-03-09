## Simple code to run SCENIC on your single cell data

The fundamental big plus of single cell seq experiments is exploring the heterogeneity of a biological sample. 
This is particularly interesting (I think) if you are focusing on a single tissue type as it undergoes change
This might be cells in a developing tissue, in culture or cells in a disease process (maybe compared to normal tissue)

When you are doing the above, scRNAsq is essentially showing you the various _transcriptional states_ of a tissue. 
How many states exist within your sample is a tricky question, and kind of centres around some sort of combination
of how diverse the sample actually is, combined with your clustering resolution and how confident you are that each
cluster you annotate as a separate state/type is actually that. 

Once you've identified these changing cell states, it makes sense to ask what is driving or regulating the changes.
In other words, **Which transcription factors are pulling the strings?**

TO keep it SCimple, SCENIC does two things:

**1)** It finds out which genes are co-expressed with TFs in the data, and writes that down in a table

**2)** It scrutinises the results to try and reduce noise (chance co-expression)

The latter step is more complex, but basically it uses databases (you need to download these and point SCENIC in their direction)
to check if the appropriate **transcription factor binding motifs** are significantly over-represented near the **transcription start site (TSS)**
of the genes picked out in step 1). This leaves you with a set of TFs and their targets, known together as a **'regulon'**

Whilst the results are still going to be noisy, it remains a powerful tool. However, as datasets get bigger and bigger, the computational demands of running
SCENIC in R are becoming overwhelming. 14k cells run on 100gb cluster takes about a week.. So they have a faster python implementation (not great if you fear snakes)

The other issue is if like me you're not really a computer person, the vignettes aren't great and contain several typos and outdated commands.

I've put this code together to solve the above problems- it allows you to run only the intensive part in python (for speed) and everything else including plotting in R.

I'm assuming an annotated Seurat object as a starting point, and use the human database (see https://resources.aertslab.org/cistarget/ for others)

```markdown
# First, install all the packages you're going to need (I'll assume you have the more amainstream ones)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
BiocManager::install(c("doMC", "doRNG"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC") 

# load libraries
library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table)
library(grid)
library(ComplexHeatmap)

# load seurat object
load("<path-to-your-object>")

# optional- subset your object to make it more manageable (eg 10k cells)
seurat = subset(seurat, cells = sample(Cells(seurat), 10000))
  
# At this point it's probably best to create a working directory to do all this work in, so do this next
dir.create("SCENIC")
setwd("SCENIC")

# download the database to this directory so that SCENIC can use it later on (they're quite large)
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
"https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) 
}
           
# Create expression matrix
exprMat <- seurat[["RNA"]]@data
exprMat<- as.matrix(exprMat)

# Create CellInfo object and save in a new directory 'int' to allow analysis by celltype later on (messy code to avoid an error later on)
cellInfo <- data.frame(seuratCluster=Idents(seurat))
cellTypeColumn <- "seuratCluster"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$CellType))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
          
# save colour info for each of your cell states/types in the above CellInfo - use the same colours as your seurat 
dimplot to make for better figures later on
colVars <- list(CellType=c("cellA"="#1F78C8", "CellB"="#ff0000", "CellC"="#33a02c", "CellD"="#6A33C2", "CellF"="#ff7f00", "CellG"="#565656"        
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")

# quickly check the colours match the UMAP
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

# now actually create your scenic object. The object class is its own thing - scenicOptions - and as you perform functions
it will write **a lot** of new files to the "/int" directory and a "/output" directory it will create
#first tell it the directory where the databases are and that they are human                    
scenicOptions <- initializeScenic(org="hgnc", dbDir="scenic", nCores=10, dbs=defaultDbNames[["hgnc"]])

# now the cell info
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds" 
                                           
# now tell it where to put the loom object it will create (eventually) from the anaylsis, and what to call that loom object
scenicOptions@fileNames$output["loomFile",] <- "output/yourloom.loom"
  
# processing the raw matrix
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)   

# save the scenicOptions object
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
```

At this stage we are ready to test for TF / gene co-expression

This can be done in R, using GENIE3 (gold standard, but takes a loooooooong time) or in python, using GRNBoost (around 10x faster (so still not quick) and essentially the same performance)

### The Python route
```markdown
# In R, this command writes the approapriate txt files for use in python in the 'int' directory from before
exportsForArboreto(exprMat, scenicOptions, dir = "int")
```

Then switch to python!
```markdown
#you may need to install some packages first. This is best done in a new environment, unless you're using a jupyter notebook (recommended for ease)
# create an environment called 'scenic'
$ conda create --name scenic-env

# activate it
$ source activate scenic-env

# now install arboreto (this gives us access to grnboost)
(arboreto-env) $ pip install arboreto

#import
import os
import pandas as pd

from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

#load in the .txt files created in R
tf_names = load_tf_names("<path-to-txt-file>")
ex_matrix = pd.read_csv("<path-to-text-file>, sep='\t')
  
#run GRNboost
network = grnboost2(expression_data=ex_matrix, tf_names=tf_names)
  
#write the output for use as the linklist in downstream SCENIC R steps
network.to_csv('linklist.tsv', sep='\t', header=False, index=False)
  ```

### The R route - slow and steady..
```markdown
# this step will take many days with modern datasets. run it in terminal to reduce the risk of crashing adn lsing everything.
runGenie3(exprMat_filtered_log, scenicOptions)

# a clever failsafe is to split the exprmat into a list, splitting the task into small chunks which are saved separately then combined. 
# Requires GENIE3 direct rather than in SCENIC- note altered command (source: https://github.com/aertslab/GENIE3/issues/1)
# the "GENIE3_linkList.RData" is the same end output as the above command in SCENIC
library(GENIE3)
genesSplit <- split(sort(rownames(exprMatrix_filtered_log)), 1:10)
lenghts(genesSplit)

for(i in 1:length(genesSplit))
{
  print(i)
  set.seed(93827)
  weightMatrix <- GENIE3(exprMatrix_filtered_log, targets=genesSplit[[i]])
  save(weightMatrix, file=paste0("GENIE3_weightMatrix_",i,".RData"))
}

# Merge results:
library(GENIE3)
linkList_list <- list()
for(i in 1:10)
{
  load(paste0("int/1.3_GENIE3_weightMatrix_",i,".RData"))
  linkList_list[[i]] <- getLinkList(weightMatrix)
}
length(linkList_list)
sapply(linkList_list, nrow)

linkList <- do.call(rbind, linkList_list)
colnames(linkList) <- c("TF", "Target", "weight")
linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
linkList <- linkList[which(linkList[,"weight"]>0),]
save(linkList, file="GENIE3_linkList.RData")
```
### Whichever route you chose, you should now have a 'linklist' in your directory, +/- individual GENIE3 outputs...
![image](https://user-images.githubusercontent.com/77628512/110475406-bef5c080-80d8-11eb-89a3-abe3a47416a6.png)

### Now we need to do the final steps back in R using this linklist as the starting point
```markdown
scenicOptions <- readRDS("int/scenicOptions.Rds")

# Build and score the GRN- step 2 mentioned at the start. worth reducing the number of cores as for unknown reason (at time of writing) this causes errors if too high
# https://github.com/aertslab/SCENIC/issues/38 https://github.com/aertslab/AUCell/issues/3 
scenicOptions@settings$nCores <- 2

# In the first step, you need to specify your link list if you used the python step
runSCENIC_1_coexNetwork2modules(scenicOptions, linkList = "<path-to-linklist.tsv>"))

# now run the rest
runSCENIC_2_createRegulons(scenicOptions)
exprMat_log <- log2(exprMat+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# write the results into a loom file (specified when you initialised the scenicOptions) to explore the results. Two weird lines added here due a bug in the code
fileName <- getOutName(scenicOptions, "loomFile")
scenicOptions@inputDatasetInfo$cellInfo = cellInfo
export2scope(scenicOptions,exprMat)
```
After this stage, everything you need for analysis is done and can be loaded at any time:

```markdown
scenicOptions <- readRDS("int/scenicOptions.Rds")
cellInfo<- readRDS("int/cellInfo.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulonsAuc(loom)
regulonsAucThresholds <- get_regulonThresholds(loom)
embeddings <- get_embeddings(loom)
```

