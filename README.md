# Source code for *Drivers of multifaceted beta diversity change in invaded stream fish communities*
## Contact information and citation

```bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```
Cite as:
> CITE


## Analysis work flow
Here, we provide code and a general workflow for calculating null model 
standardized alpha diversity, beta diversity, and local contribution to beta 
diversity (LCBD) values. Standardized effect sizes (SES) were calculated for 
two species pools (contemporary and native only), as well as for the change in 
diversity between species pools. See the manuscript for more detailed
methodology and justification. This workflow is designed to run using R on 
local machines with the more expensive calculations ran using high performance
computing clusters via the slurm interface and shell scripts. R and shell 
scripts are numbered in order of workflow.  Scripts 1-14 are coded to be general 
so that users can calculate SES using their own diversity values. We used this 
framework to create the diversity values used the analyses in this manuscript. 
Scripts 15-17 are less general and are provided with the goal of transparency 
and result replication. We provided both raw and formatted null model data, 
explanatory variables, and the sources for publicly available explanatory data 
to replicate analyses from Script 11 onward. We are unable to directly provide 
the raw community or trait data that facilitate replication of scripts 1-10, but 
provide additional code and workflow information regarding the filtering and 
formatting of the diversity input data (i.e., community, trait, phylogeny). The 
final phylogenetic tree, and information for accessing community and trait data 
can be found [here](#diversity-input-data).


### Required Software
**R version**: 4.5.0

R packages for null model workflow
* ```'ade4'``` version: 1.7.23
* ```'adespatial'``` version: 0.3.28
* ```'ape'``` version: 5.8.1
* ```'BAT'``` version: 2.11.0
* ```'DescTools'``` version: 0.99.60
* ```'dplyr'``` version: 1.1.4
* ```'ggplot2'``` version: 4.0.2
* ```'matrixStats'``` version: 1.5.0
* ```'parallel'``` version: 4.5.0
* ```'purrr'``` version: 1.1.0
* ```'tibble'``` version: 3.2.1
* ```'VGAM'```  version: 1.1.14

R packages for mansucript analysis replication
* ```'sf'``` version: 1.0.20
* ```'StreamCatTools'``` version: 0.10.0
* ```'tidyr'``` version: 1.3.1
* ```'vegan'``` version: 2.6.10

R packages for Diveristy input prep
* ```'fishtree'``` version: 0.3.4’
* ```'picante'``` version: 1.8.2
* ```'RRphylo'``` version: 3.0.2
* ```'stringr'``` version: 1.5.1


### Create file directories
First, download the ```Scripts``` folder (all users), all data from 
[Zenodo Repository]() (users replicating results), and all explanatory variables 
(users replicating results) from data sources listed [here](#analysis-data). 
Next, users will need to create the below file directory to store diversity 
input data (e.g. community, trait, phylogeny), formatted high performance 
computation input data, and the resulting diversity outputs.

```bash
├── Scripts
├── Diversity Input Data
│   ├── mod_com_diversity_input.rds
│   ├── his_com_diversity_input.rds
│   ├── trait_diversity_input.rds
│   ├── phylo_tree.rds
├── HPC
│   ├── 00_null_model_effect_size_function.R
│   ├── 01_ses_batch.sh
│   ├── 01_ses_batch.R null_model_effect_size_function.R
│   │── null_out*
│   │── obs_out*
│   │── ses_inputs
│   └── ses_outputs*
├── Diversity Output Data
├── Analysis_data
└── Results
```
Place community data for both species pools (contemporary - "mod" and native - "his"), trait data, and phylongetic trees in ```Diversity Input Data``` with the following names: 
* ```mod_com_diversity_input.rds```: Community data for contemporary species pool. Dataframe or matrix with rows for sites and columns for species. Can contain an optional column ```HUC_12``` which represent regions and can be used to define regional species pools.
* ```his_com_diversity_input.rds```: Community data for native-only species pool. Same structure as ```mod_com_diversity_input.rds```
* ```trait_diversity_input.rds```: Trait data for all species in community data. Dataframe or matrix with rows for species and columns for traits.
* ```phylo_tree.rds```: Phylogenetic tree for all species in community data

<ins>NOTE:</ins> We cannot provide raw community or trait data used in the manuscript without completed data requests, but we do provide the phylogenetic tree and scripts used to format and filter the trait and community data. See [Diversity Input Data](#diversity-input-data) for more information

<ins>NOTE:</ins> If replicating the results of analysis (Script 11 onward) by downloading data from the [Zenodo respository](), the .zip files will create duplicate folders. It is recommended to download all data before creating new file structures.


### Prepare input data - perform on local machine
To best take advantage of high performance computation, we need to format our data in a manner that allows for parallel processing of diversity metrics. As we need to estimate 999 null iterations of each diversity metric and species pool (native only and contemporary), we want to prepare input data that can run iterations simultaneously. For beta diversity null models, we need to generate 999 shuffled trait matrices and phylogenetic trees. For native alpha diversity, we need to generate 999 random community matrices. A benefit to HPC is that processes can be ran on a high number of cores among multiple computer nodes. To take advantage of HPC, we need to divide the list of null traits, trees, and/or communities into chunks based on CPU and memory limits per each node. These chunks can be ran on separate computer nodes. This reduces memory requirements within nodes and allows for better queue times. 

The following script imports ```mod_com_diversity_input.rds```, ```his_com_diversity_input.rds```, ```trait_diversity_input.rds```, and ```phylo_tree.rds```. These data are then used to create reduced functional space using PCoA and trims trees to fit community data. Both observed and null input data lists for alpha and beta diversity are then created and divided into chunks to run on separate HPC nodes.
* ```01_null_input_creation.R```

After running the above script, upload the entire ```HPC_data``` directory to the high performance cluster storage.

<ins>NOTE:</ins> This script only contains code to shuffle communities using taxa-swap and a regionally-constrained taxa-swap null model algorithms, which are ran using functions called from ```null_model_algorithms.R```. While the algorithms provided were best suited for functional and phylogenetic beta diversity, they may not be suited for all null model purposes. Users should research the best model for their usage and modify ```01_null_input_creation.R``` and ```null_model_algorithms.R``` accordingly.


### Calculate observed diversity values - perform using HPC
The below shell scripts call in their corresponding R scripts to estimate observed beta diversity for the contemporary and native only species pools using one high performance computer nodes for each time step. Alpha diversity scripts estimate alpha diversity of only native species using one computer node.
* ```02_obs_tax_beta.sh``` - ```02_obs_tax_beta.R```
* ```03_obs_fun_beta.sh``` - ```03_obs_fun_beta.R```
* ```04_obs_phy_beta.sh``` -```04_obs_phy_beta.R```
*	```05_obs_fun_alpha.sh``` - ```05_obs_fun_alpha.R```
*	```06_obs_phy_alpha.sh``` - ```06_obs_phy_alpha.R```


### Calculate null diversity values  - perform using HPC
The below shell scripts call in their corresponding R scripts to estimate null iterations of beta diversity for the contemporary and native only species pools. Due to the high memory usage of kernel density functional beta diversity metrics, we were only able to estimate 2 iterations (2 cores) and 18gb of ram per computer node for a total of 2000 high performance nodes. For phylogenetic beta diversity, we were able to estimate 37 iterations simultaneous per computer node using 45gb of ram per node for a total of 6 high performance nodes. For alpha diversity metrics, we estimated only null iterations for the native species pool and this required less than half the resources of beta diversity null models. 

*	```07_null_fun_beta.sh``` - ```07_null_fun_beta.R```
*	```08_null_phy_beta.sh``` - ```08_null_phy_beta.R```
*	```09_null_fun_alpha.sh``` - ```09_null_fun_alpha.R```
*   ```10_null_phy_alpha.sh``` - ```10_null_phy_alpha.R```

<ins>TIP:</ins> Number of nodes, cores per node, and memory per node will vary based on number of sites and methodology. Shell scripts can be edited to adjust these settings accordingly. For example: 1000 nodes: ```#SBATCH --array=1-1000```, 2 CPU cores per node: ```SBATCH --cpus-per-task=2 ```, 18gb ram per node :```#SBATCH --mem=18gb```. We recommend that users experiment with memory and CPU requirements with smaller number of null iterations before running full job.


### Calculate effect sizes - perform using HPC
Observed and null model outputs were exported in multiple files, reflecting the multi-nodal processing. We need to consolidate these files into single files that contain a list of each null iteration.  The below shell scripts run their respective R script to consolidate null model outputs into single files, separate native and contemporary species pool values, estimate the difference in diversity between contemporary and native pools (delta), and export a single file for each diversity metric (i.e., alpha, total beta, replacement, richness difference, LCBD) for the contemporary species pool, native only species pool, and delta values.

* ```11_beta_null_model_prep.sh``` - ```11_beta_null_model_prep.R```
* ```12_alpha_null_model_prep.sh``` - ```12_alpha_null_model_prep.R```

The lists from the functions above contain the observed values and a list of null iterations for each metric and species pool. The next shell script will run the respective R script, which estimates standardize effect sizes (SES) of each single metric. This step involves calling in a SES function created for this project: ```null_model_effect_size_function.R```. This function is flexible and takes a range of input formats such as dataframes, vectors, matrices, and distance objects, and maintains this format in the exported values. The function estimates standardize effect sizes in the traditional z score method (SES). Additionally empirical p-values, and p-value based effect sizes (ES) are calculated. Finally the function reports optional diagnostic metric to assess if null distributions are symmetrical and normal. asymmetrical null distributions should be assessed using empirical p-value based effect sizes rather than z-score based SES. See [Botta-Dukát (2018](https://doi.org/10.1556/168.2018.19.1.8) for more information on selecting SES or p-value based ES.

* ```13_batch_ses.sh``` - ```13_batch_ses.R```

After running the above scripts, download the entire ```HPC_data``` directory to local machine.


### Summarize null model results - perform on local machine
Batch SES processing results in a single file for each diversity metric. The following R script compiles and formats the resulting SES, ES, and diagnostic stats across files. Here, we can visualize the normality diagnostics and determine if we need to use SES or ES values for further analysis. The script also visualizes beta diversity and LCBD change values and creates plots for Figure 4 in the manuscript. Finally, the script prepares and exports native alpha/LCBD and delta LCBD for use in the manuscripts analyses.

* ```14_ses_comp.R```

<ins>NOTE:</ins> This is final step for creating null model effect sizes for alpha, beta, and LCBD diversity values. All remaining workflow is for replication of the manuscripts results. To facilitate replication, we provide the raw [observed](#observed-diversity-data), raw [null iterations](#null-iterations), [summarized null model results](#summarized-null-model-outputs), and the [formatted alpha and LCBD values](#diversity-output-data) used for manuscripts analyses.


### Redundancy analaysis - perform on local machine
The below script prepares the explanatory variables for used in redundancy analysis (RDA) of multidimensional changes in LCBD. The script loads in data for nonnative origin-based invadedness, propagule pressure, abiotic habitat characteristics, habitat alteration, and native alpha diversity and LCBD. This script also summarizes community invadedness by species origin and creates the table for Appendix 7. We make all data for explanatory variable available or explain how to access the data. See [Analysis Data](#analysis-data) for more information.
* ```15_rda_predictor_prep.R```

The below script conducts forward selection, redundancy analysis, and variance partition for both raw (observed) and null model standardized (ES) change in LCBD values. Script also contains code to export tables and plots for figure 6 and appendix 8 of the manuscript.
* ```16_rda_varpart.R```


### Spatial plotting - perform on local machine

The following script aggregates the delta LCBD values from stream segment to HUC6 basin resolution to visualize spatial patterns of changes in LCBD of total beta diversity. Shapefiles are exported into QGIS for mapping and formatting to create Figure 5 in manuscript.
* ```17_lcbd_spatial.R```


### Helper functions
* ```null_model_algorithms.R```: Contains algorithms to randomize community, trait, or phylogenetic data for null model analysis
* ```diversity_batch_functions.R```: Contains functions that estimate multiple iterations of diversity metrics using parallel computation.
* ```null_model_effect_size_function.R```: Function to summarize null distributions and estimate effect sizes. Function estimates standardize effect sizes in the traditional z score method (SES), empirical p-values, and p-value based effect sizes (ES). Finally the function reports optional diagnostic metric to assess if null distributions are symmetrical and normal.

## Diversity Input Data
We are not able to publicly provide all the data necessary to replicate the multidimensional diversity data. The fish occurrence data used to estimate multidimensional diversity metrics were obtained through data sharing agreements with United States governmental agencies. While these raw data are not directly available from the authors for redistribution due to data sharing agreements, they can be accessed through formal requests to the agencies listed in Appendix 1. Individuals with completed data requests may contact the corresponding author for harmonized versions of the data. Phylogenetic data was obtained from a variety of publicly available datasets with citations found in Appendix 3. We provided a harmonized phylogenetic tree that contains all species in our community dataset. Trait data were obtained from public and private datasets and harmonized data cannot be shared without permission. Despite not being able to share all data, we provide R scripts used to harmonize the community, phylogenetic, and trait data used to estimate the multidimensional alpha and beta diversity values, which we do make publicly available.


### Community  data
Harmonized community data is avalable upon request after written approval from each data source in appendix 1 of manuscript. Community data should be stored: ```Diversity Input Data/```. While we cannot publicly share raw community data, we provide the list of species in the analysis for use in the trait and phylogenetic data preparation scripts.
* ```full_species_list.csv```: List of all species in the community data set before filtering, used for construction of phylogenetic super tree and compilation of trait data.
* ```filtered_species_list.csv```: List of all species in the community data set after filtering, used for to trim tree and compile trait data. 
* ```community_data_prep.R```: Compiles stream fish community data at the stream segment scale for the conterminous United States (US) and to filter data to only include comparable surveys based on sampling methods. This code also rarefies stream segments based on hydrological region to ensure consistent sampling densities across the study extent. Finally, stream segments are split into two species pools, a contemporary that includes all species records, and native only species pool, which represents communities before nonnative introductions. <ins>NOTE:</ins> Running this script requires data not publicly available without data requests.

### Phylogenetic data
* ```phylo_tree.rds```: Phylogenetic super tree created using the methods of [Castiglione et al.,  (2022)](https://doi.org/10.1111/pala.12588) and ```phylo_data_prep.R```. This tree includes all 815 species from the community dataset used in the analysis. This tree used the [Fish Tree of Life (Rabosky et al., 2018)](https://fishtreeoflife.org/) as a backbone (Contains ~89% of species) with phylogenetic information supplemented by [TimeTree (Kumar et al., 2022)](http://www.timetree.org/), and several genus specific trees found in Appendix 3. <ins>NOTE:</ins> Phylogenetic super trees are constructed from trees created with differing methodology and caution should be used when making evolutionary inferences with super trees. However they are adequate for capturing general phylogenetic relationships among species.
* ```phylo_data_prep.R```: Compiles phylogenetic data from multiple trees into a super tree that includes all 449 fish species in the analysis. This script loads in a compressive fish phylogeny from the [Fish Tree of Life](https://fishtreeoflife.org/), and uses downloaded data from [TimeTree](http://www.timetree.org/) and other publicly available datasets to fill in missing species. The output of this script was used to estimate missing trait values and to calculate phylogenetic beta diversity. <ins>NOTE:</ins> Running this script requires users to download all phylogenetic trees from their original data sources found in Appendix 3 of manuscript.

### Trait data
Raw trait data obtained from publicly available datasets and upon request should be stored ```Diversity Input Data/```.
* ```trait_data_prep.R```: Compile trait data for all 449 fish. This script uses trait data from two trait databases and fills in missing trait values based on phylogenetic relationships and literature review. The output of this script was used to calculate functional beta diversity. <ins>NOTE:</ins> Running this script requires users to download the publicly available data set from [Frimpong & Angermeier (2009)](https://www.sciencebase.gov/catalog/item/5a7c6e8ce4b00f54eb2318c0) and request the database from [Giam & Olden (2016)](https://doi.org/10.1111/geb.12475).

## Observed diversity data
Observed alpha and beta diversity data can be downloaded from ```obs_out.zip``` at the [Zenodo repository]() and unzipped into ```HPC_data/``` directory. These data 
used to estimate null model empirical effect sizes, used in analyses, and used for summary statistics.

### Observed beta diversity data
```his_fun_beta_obs.rds``` ```his_phy_beta_obs.rds``` ```his_tax_beta_obs.rds``` ```mod_fun_beta_obs.rds``` ```mod_phy_beta_obs.rds```
```mod_tax_beta_obs.rds```

Files contain lists of pairwise beta diversity distance objects and local contribution to beta diversity (LCBD) data frames. Files contain values for the contemporary (mod) and native (his) species pools, and taxonomic (tax) functional (fun), and phylogenetic (phy) diversity facets. All files have the following structure:
* ```$beta```: list of 3 containing Sorensen pairwise beta diversity distance objects
    * ```$Btotal```: distance object of total beta diversity
    * ```$Brepl```: replacement beta diversity component
    *  ```$Brich```: richness difference beta diversity component
*  ```$LCBD```: Data frame of 3 columns containing LCBD values. Column names are specific to diveristy facets, and these differences are denoted with *X*:
    * ```X_Btotal```: LCBD of total beta diversity in facet *X*
    * ```X_Brepl```: LCBD of replacement component in facet *X*
    * ```X_Brich```: LCBD of richness difference component in facet *X*

### Observed alpha diversity data

```his_fun_alpha_obs.rds``` ```his_phy_alpha_obs```

Files contain named numeric objects with observed functional richness (volume of kernel density hypervolume) or phylogenetic richness (number of branches in phylogenetic tree). Names refer to COMID identifiers for each stream segment via the [National Hydrography Dataset Plus version 2](https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data#Download)

## Null iterations
The raw null model iteration data for beta and alpha diversity can be downloaded from ```null_out.zip``` at the [Zenodo repository]() and unzipped into the ```HPC_data/``` directory. These data are the result of the randomization of traits and phylogenies (beta) or randomization of communities (alpha). These data will be compiled and used to create null distributions used to create null model standardized diversity values.

Null model iterations are divided into separate files for alpha and beta diversity based on diversity facet. Within diversity facets, iterations are broken into chucks based on the number of HPC nodes were used to create the data. The null model files are named in the following format
> *pool*\_*facet*\_*metric*\_null\_*iteration*.rds
* *pool*: species pool (contemporary = mod; native = his)
* *facet*: diversity facet (taxonomic = tax; functional = fun; phylogenetic = phy)
* *metric*: alpha or beta diversity
* *iteration*: range of iterations included in file. e.g., 001-002

All files contain a list with items for each null iteration. Each beta diversity iteration contains a list with the same structure as [observed beta diversity](#observed-beta-diversity-data) and each alpha diversity iteration contains a named numeric object with the same structure as [observed alpha diversity](#observed-alpha-diversity-data).

   
## Summarized null model outputs
The summarized null model analyses results for all diversity metrics can be downloaded from ```ses_out.zip``` at the [Zenodo repository]() and unzipped into the ```HPC_data/``` directory. These data are used to evaluate the properties of the null distributions to choose between standardized effect sizes and empirical effect sizes, and can be merged together for plotting and to create the final dataframes for the use in analyses.

Each diversity metric has its own file name in the following format:
> *facet*\_*pool*\_*metric*\_ses_out.rds
* *facet*: diversity facet (taxonomic = tax; functional = fun; phylogenetic = phy)
* *pool*: species pool (contemporary = mod; native = his)
* *metric*: diversity metric (Btotal = total beta diversity; Brepl = replacement component; Bric = richness difference component; local contribution to beta diversity = LCBD; alpha = alpha diversity)

All files contain a list with the following structure:
* ```$obs```: observed diversity values
* ```$null_mean```: means of null distributions
* ```$null_sd```: standards deviation of null distributions
* ```$ses```: z-score based standardized effect sizes: ```(obs - null_mean)/null_sd)```
* ```$empirical_pvalue```: proportion of null distribution with higher values than observed value
* ```$empirical_es```: Empirical p-value based effect sizes ```probit(1-empirical_pvalue)```
* ```$skew```: skewness of null distribution
* ```$kurt```: kurtosis of null distribution

The structure of the values in each list item are dependent on diversity metric: 
* beta diversity components (Btotal, Brepl, and Brich): distance objects,
* LCBD: data frame with rows for each site and a column for LCBD of Btotal, Brepl, and Bric 
* alpha: named numeric objects.


## Diversity Output Data
The combined observed and summarized effect size data for native alpha, native LCBD, and delta LCBD can be downloaded from ```Diversity Output Data.zip``` at the [Zenodo repository]() and unzipped into the working directory. These data can be used to replicate the spatial, redundancy, and variance partitioning analyses.

* ```delta_lcbd.rds```: Observed and null model empirical p-value based effect sizes (ES) of taxonomic (no ES), functional and phylogenetic changes in local contributions to beta diversity (LCBD) between the contemporary and native species pools. These data are used as the response variables for the manuscript's main analyses. Data frame consisting of the following columns:
   * ```COMID```: Unique identifier for each stream segment via the [National Hydrography Dataset Plus version 2](https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data#Download)
    * ```fun_Btotal```: Change in functional LCBD of total beta diversity
    * ```fun_Brepl```: Change in functional LCBD of the replacement component
    * ```fun_Brich```: Change in functional LCBD of the richness difference component
    * ```phy_Btotal```: Change in phylogenetic LCBD of total beta diversity
    * ```phy_Brepl```: Change in phylogenetic LCBD of the replacement component
    * ```phy_Brich```: Change in phylogenetic LCBD of the richness difference component
    * ```tax_Btotal```: Change in taxonomic LCBD of total beta diversity
    * ```tax_Brepl```: Change in taxonomic LCBD of the replacement component
    * ```tax_Brich```: Change in taxonomic LCBD of the richness difference component
    * ```fun_Btotal_es```: Empirical effect size of change in ```fun_Btotal```
    * ```fun_Brepl_es```: Empirical effect size of change in ```fun_Brepl```
    * ```fun_Brich_es```: Empirical effect size of change in ```fun_Btotal```
    * ```phy_Btotal_es```: Empirical effect size of change in ```phy_Btotal```
    * ```phy_Brepl_es```: Empirical effect size of change in ```phy_Brepl```
    * ```phy_Brich_es```: Empirical effect size of change in ```phy_Brich```
* ```native_lcbd.rds```: Observed and ES values of taxonomic (no ES), functional, and phylogenetic local contributions to beta diversity (LCBD) for the native species pool. These data are used as explanatory variables in main analyses. Data frame consisting of the same structure as ```delta_lcbd.rds```, but LCBD values are for the native species, not a change over time.
* ```native_alpha.rds```: Observed and ES values of functional and phylogenetic richness for the native species pool. These data are used as explanatory variables in main analyses. Data frame consisting of the following columns:
    * ```COMID```: Unique identifier for each stream segment.
    * ```fun_his_alpha```: Functional richness measured as the volume of kernel density hypervolumes.
    * ```phy_his_alpha```: Phylogenetic richness measured as the number of phylogenetic tree branches.
    * ```fun_his_alpha_es```: Empirical effect size of ```fun_his_alpha```
    * ```phy_his_alpha_es```: Empirical effect size of ```phy_his_alpha```

## Analysis data
We directly provided the data or the sources to the data required to replicate the redundancy analysis, variance partitioning, and spatial analysis. The response variables ```delta_lcbd.rds``` and biotic explanatory variables ```native_lcbd.rds``` ```native_alpha.rds``` for analyses can be found in [Diversity Output Data](#diversity-output-data). Information on how to access the remainder of the explanatory variables is listed below.

### Origin-based invadedness
We estimated the effects community invadedness by nonnative species of different geographic origins on changes in LCBD. We used the raw community data estimated total species richness, native species richness, and the richness of three classes of nonnative species using origin-base definitions introduced by [Thompson et al., (2025)](https://doi.org/10.1111/geb.13951). We are not able to publicly share the raw community data but provided the summarized richness values and r script ```origin_invaded_prep.R``` used to generate the data. The R script can be found in```analysis_data/``` and data can be downloaded from ```origin_invaded.rds``` at the [Zenodo repository]() and unzipped into ```analysis_data/```.
* ```origin_invaded.rds```: Species richness based on origin based native status, community invadedness, and region to site bridge data. Data frame with the following columns:
    * ```HUC_12```: Unique identifier for each subwatershed (HUC_12) via the [National Watershed Boundary dataset](https://www.usgs.gov/national-hydrography/watershed-boundary-dataset). Can be used to link stream segment to any level of watershed organization (i.e., HUC2-12).
    * ```COMID```: Unique identifier for each stream segment.
    * ```nTaxA```: Number of species native to HUC8
    * ```prov```: Number of provincially nonnative species. These are species nonnative to the province (HUC8), but native to the region (HUC2).
    * ```reg```:  Number of regionally nonnative species. These are species nonnative to region, but native to the realm (continent).
    * ```extr```:  Number of extra-realm nonnative species. These are species nonnative to the realm.
    * ```tot_sp```: Total number of species (native and nonnative).
    * ```invProv```: Community invadedness (i.e., Proportion of total species richness) by provincially nonnative species.
    * ```invReg```: Community invadedness by regionally nonnative species.
    * ```invExtr```: Community invadedness by extra-realm nonnative species.

### Propagule pressure
We estimated the effects of recreational fishing demand on changes in LCBD using a metric created by [Mazzotta et al. (2015)](https://doi.org/10.1016/j.ecolecon.2015.09.018) and [Davis and Darling (2017)](https://doi.org/10.1111/ddi.12557). These data can be downloaded from the [EPA EnviroAtlas](https://www.epa.gov/enviroatlas) and placed into the ```analysis_data/``` directory.

### Abiotic habitat characteristics 
We estimated the effects of elevation, 30-year mean temperature, base flow index, and upstream watershed area. These data are all downloaded in ```15_rda_predictor_prep.R``` using the ```'StreamCatTools'``` package. See [here](https://github.com/USEPA/StreamCatTools) for more information on the ```'StreamCatTools'``` package.

### Abiotic habitat alteration
We estimated the effects of hydrological alteration on changes in LCBD using an index of hydrological alteration developed by [McManamay et al.,(2020)](https://doi.org/10.1038/s41597-022-01566-1). These data can be downloaded [here](https://zenodo.org/record/5839011) and placed into the ```analysis_data/``` directory.

### Biotic Factors
We estimated the effects of native alpha diversity and LCBD on changes in LCBD. Information on obtaining these values can be found in [Diversity Output Data](#diversity-output-data) in the following files: ```native_lcbd.rds``` and ```native_alpha.rds```.

### Spatial data
We mapped the spatial distribution of changes in LCBD by aggregating data to the HUC6 basin level. Doing so requires ```origin_invaded.rds``` to bridge between COMIDs and HUC6. Once bridged, LCBD values can be attached to HUC6 shapefiles downloaded from the [USGS National Hydrography Watershed Boundary Dataset](https://www.usgs.gov/national-hydrography/watershed-boundary-dataset). Place downloaded shapefiles in the ```analysis_data/``` directory.
# Beta_diversity_change_drivers
