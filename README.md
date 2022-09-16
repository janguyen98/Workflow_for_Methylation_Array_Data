# Rheumatoid arthritis methylation microarray data analysis
General workflow for analysis of methylation microarray data of rheumatoid arthritis samples.

## General Workflow

### 1. Loading the data

R packages used:

* minfi
* limma
* IlluminaHumanMethylation450kanno.ilmn12.hg19
* IlluminaHumanMethylation450kmanifest
* missMethyl
* minfiData
* DMRcate
* MethylToSNP
* FlowSorted.Blood.450k
* sva
* glmnet

Metadata is loaded in as an R object, then written out to CSV. `read.metharray.sheet` reads in metadata CSV as a table which is then used to create an `RGChannelSet` object for further analysis.

### 2. Quality control of samples

Detection p-values are calculated for each probe. Samples whose detection p-values are higher than 0.05 for at least 10% of the probes are dropped.

### 3. Proprocessing

Typically, functional normalization applied and is a common preprocessing method used for global methylation differences, such as in cancer studies or in studies with different tissue types. Other preprocessing methods include:

* `preprocessRaw`: no normalization applied
* `preprocessQuant`: uses stratified quantile normalization when cell types are the same; forces distributions of the populations to be the same, removing all variation
* `preprocessFunnorm`: functional normalization extends the idea of quantile normalization by removing unwanted variation by regressing out variability explained by the control probes present on the array; applicable for cases where global changes are expected such as in cancer-normal comparisons or tissue differences.
* `preprocessNoob`:  NOOB normalization applies a background correction based on the negative control probes
* `preprocessSWAN`: normalizes Type I probes to Type II probes

#### 3.1 Batch correction

Common unsupervised batch correction methods are:

* **Surrogate variable analysis (SVA)**: estimates hidden covariates using surrogate variables; may remove unknown yet important biological differences between samples
* **ComBat**: adjusts for batch effects when batches are known; can introduce false signal if study design is unbalanced
* **Removing Unwanted Variation (RUV)**: uses factor analysis on control genes; must select the control genes
* **Hidden Covariates with Prior (HCP)**: batch effects (both known and unknown) are modelled as a linear combination of known covariates
* **Batch Effect Signature Correction (BESC)**: uses a reference dataset to learn variation due to unknown technical differences and apply the correction to other datasets

### 4. Data exploration

MDS plots are used to examine how samples are clustered. Covariates are then appropriately adjusted for during probe filtering or during linear model fitting.

### 5. Probe filtering

Probes are flagged or filtered out to remove poor-quality probes or those that will skew results of the analysis. Those that are filtered out include:

* probes with detection p-values more than 0.05 in at least 20% of samples
* probes on sex chromosomes
* cross-hybridizing probes (as identified [here](https://github.com/sirselim/illumina450k_filtering))
* probes with SNPs (use *MethyltoSNP* package or *minfi's* `dropLociWithSnps`)

### 6. Differential methylation calculation

Beta values and M values are calculated using `minfi`. Beta values are easier to interpret biologically, since values range from 0 to 1 (0 being the site is unmethylated and 1 being the site is methylated). M-values are logit-transformed beta values for use in statistical analysis ([Du et. al. 2010](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587#article-info)).

Analyses typically identify differentially methylated positions (DMPs), differentially methylation regions (DMRs), or differentially variable positions (DVPs). If using `limma`, a model matrix is constructed to include all covariates, then a linear model is fitted to find DMPs. 

### 7. Functional annotation and gene set enrichment analysis (GSEA)

After identifying statistically significant positions or regions exhibiting differential methylation, functional annotation is performed to identify any important pathways relevant to the analysis:

* DAVID bioinformatics
* Gene Set Enrichment Analysis (GSEA)
* g:Profiler

### 8. Comparison with Methyl-Seq data

Pairing methylation microarray data with methylation sequencing data could illuminate new differentially methylated sites not captured by the microarray.

### 9. Classification model construction

A classification model is built to help identify a list of candidate biomarkers. Types of models include lasso regression, elastic net model, and ridge regression. Read about these types of models [here](https://www.datacamp.com/tutorial/tutorial-ridge-lasso-elastic-net). Lasso or elastic net models are preferred, both of which perform feature-reduction to select only the most important probes used in classification.

### 10. Candidate biomarker selection and methylation risk score (MRS)

Candidate biomarkers are chosen following classification model construction, and methylation risk scores are calculated by multiplying the coefficiency of each biomarker by its methylation value, then taking the sum of the product. MRS is standardized by subtracting the mean MRS from each MRS, then dividing the sd MRS. Calculation of standardized MRS differs when calculating for the training set vs. validation set.

### 11. MRS validation and threshold

An MRS above a certain threshold confers high risk of RA vs. average risk of RA.
