# NetID
The pipeline for benchmarking NetID with different imputation methods on any UMI-based scRNA-seq datasets

## Setting and Usage
---------------------------
### Step1. Set the docker envir for benchmarking
The NetID pipeline follow the benchmark framework BEELINE, the BEELINE pipeline interfaces with the implementations of various algorithms through Docker containers, users could check the tutorial https://murali-group.github.io/Beeline/ for details

### Step2. Loading and Processing Datasets
The NetID require input scRNA-seq dataset is UMI-based as it would model the gene expression distribution follow negative-binomial (NB) distribution, to this end, the NetID would check if each gene of input expression matrix follows NB distribution. Running the following command in terminal to load the dataset

```
Rscript /Rscript/NB_test.r hematopoiesis.rds /Path/to/Save/scRNA-seq/Blood
```
This script would test if your data follow the NB distribution, through calculate the concordance between NB predict zero proportion in datasets and real zero proportion

![Overview of NB](https://github.com/WWXkenmo/NetID/blob/Figure/figure_readme/NB_fit.jpg)

### Step3. Creating GEP (raw and netID)
The next step of NetID is to use RaceID to create new object, perform VarID, prunning KNN graph and make umap plots. After above steps, NetID would use SCENT to evaluate differentiation entropy of each cell types to determine root cell, and run DiffusionMap to build cell differentiation trajectory and predict pseudo-time, after that NetID would use general additive model to fit each genes with pseudo-time to find variable genes on the trajectory manifold. Finally, NetID would used geosketch to sample seed cells, and aggregate each seed cells with their neighbor cells and return final GEP. All above processing could be perform through on line code.

```
Rscript /Rscript/GRN_build_pip.r hematopoiesis.rds 30 200 TRUE 3 blood DiffusionMap 10090 /Path/to/save/VarID/ /out/put/dir/inputs
```
This script require following parameters: \
**input dataset**: raw count gene expression matrix, require .rds data type \
**The number of nearest neighbor**: VarID would create KNN-graph, and this parameter determine the K \
**The number of sampled seed cell** \
**Identity mode**: if each seed cell does not connect with each other \
**The number of minimual nn**: The minimal number of NN of seed cells \
**The name of dataset** \
**Trajectory Methods**: user could select "DiffusionMap" or "Slingshot", if user doesnt want to perform trajectory analysis, just type "NO_TI" \
**The species number**： human:9060 mouse:10090 \
**The path to save varID object** \
**The path to save resulted GEP** \
![DiffusionMap](https://github.com/WWXkenmo/NetID/blob/Figure/figure_readme/flowchart.png)

### Step4. Creating GRN and evalution
#### Gene filtering
We used `generateExpInputs.py` to filter the expression matrices according to the top500 most variable genes and significantly varying TFs (Bonferroni corrected p < 0.01). This is an example command that we applied for raw and netID expression matrices (GEP). We obtain a subset of ExpressionData.csv called 500HVG_TFs-ExpressionData.csv with the mentioned filter criteria and the corresponding STRING network filtered on the same subset of genes.(User could filter each type of network)
(Create 500HVG_TFs_raw at First)
```
python generateExpInputs.py \
-e=./inputs/ESC/raw/ExpressionData.csv \
-g=./inputs/ESC/raw/GeneOrdering.csv \
-f=./Networks/mouse/STRING-network.csv \
-i=./Network/mouse-tfs.csv \
-p=0.01 \
-c \
-n=500 \
-t \
-o=/mnt/data1/weixu/NetID/Beeline-master/inputs/blood/500HVG_TFs_raw/
```
#### Creating GRN (PIDC, GENIE3 and GRNBOOST2)
Follow the tutorial of BEELINE to create config file, we build the GRN through following command
```
python BLRunner.py --config ./config-files/blood/config_ESC_500hvg.yaml
```
However, running the GENIE3 implements in BEELINE need to cost a very long time, we recommand user to build GENIE3 network through following command
```
Rscript /Rscript/GENIE3_pip.r /Path/to/save/GEP/inputs/blood/500HVG_TFs_raw/ExpressionData.csv /Path/to/outdir/outputs/blood/500HVG_TFs_raw/
```
#### Calculation Early Precision Ratio (EPR)
The user could calculate EPR through only one command (GENIE3 example)
```
Rscript /Rscript/GRN_eval.r /Path/to/network/outputs/blood/500HVG_TFs_raw/GENIE3/rankedEdges.csv /Path/to/NN_list.csv /Path/to/filtered/network/inputs/blood/500HVG_TFs_raw/ /Path/to/return/out/outputs/blood/500HVG_TFs_raw/GENIE3/
```
![EPR](https://github.com/WWXkenmo/NetID/blob/Figure/figure_readme/blood.jpg)
