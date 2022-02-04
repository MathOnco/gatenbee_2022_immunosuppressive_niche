[![DOI](https://zenodo.org/badge/455645968.svg)](https://zenodo.org/badge/latestdoi/455645968)
# Modellng tumor-immune interactions in colorectal cancer

## Description

The code here was used to simulate tumor-immune eco-evolutionary dynamics in the origins and progression of colorectal cancer. There are two models, the first of which simulates competition between fixed immune escape stratgies. The code for this model is in no_mutation_model.R. The file mutation_model.py simulates competition between immune escape strategies in a tissue that starts off in healthy colon and includes mutation. As mutations accumulate tumor cells acquire not only immune escape mechanims, but also driver mutations that allow the cell to evolve into an adenoma cell and then a carcinoma cells. However, each mutation is acommpanied by a neoantigen that increases immune predation. There is therefore selection for tumor cells that sufficiently strong immune escape mechanisisms and low antigenicity. Please see our publication "Immunosuppressive niche engineering at the onset of human colorectal cancer" for the full details of each model.

## Installation

no_mutation_model.R is written in the R programming language. It has been tested with R version 3.6.0, and requires the following packages: ggpubr (0.2.5), ggsci (2.9), AtmRay (1.31), dplyr (0.8.5), cowplot (1.0.0, viridis (0.6.2), colorspace (1.4.1), pals (1.6).

mutation_model.py was written in Python (3.7.12), and requires pandas (0.25.1), numba (0.52.0), numpy (1.21.2), and scipy (1.7.1)

## Usage

### No mutation model

no_mutation_model.R simulates competition between fixed immune escape strategies, and was used to generate the results shown in Figure 1. A parameter grid is created to generate the parameter combinations that will compete strategies with different antigenicities and "strenghts". Here, strength refers to the amoung of immune blockade (φ) and immune suppression (σ), and antigenicity (γ) determines the immune kill rate, each of which is between 0 and 1. The number of strategy strengths to sweep over is determined by the `nstrat`, while the number of antigenicity values to sweep over is determined by the `nant` parameter. To run this model, simply run:

    Rscript no_mutation_model.R

### Mutation model

mutation_model.py simulates tumor-immune eco-evolutionary dynamics, from normal healthy tissue to adenoma to carcinoma. Similar to the no mutation model, cells have phenotypes that are combinations of immune blockade (φ), immune suppression (σ), and antigenicity (γ). However, cells aslo have driver mutations that allow them to divide more rapidly and invade. This code conducts one "run" given a set of parameters, and is intended to be used via the command line so that many runs can be conducted simultaneosly on a cluster. The command line arguments are:

* `dst_dir` (directory in which to save the results)
* `f_prefix` (file prefix for all results)
* `bp_mutation_rate` (base pair mutation rate)
* `blockade_protection` (degree to which the immune kill rate is decreased when the immune blocakde strategy is acquired)
* `immunosuppression` (amount of immunosuppression provided by immunosuppressive cells)
* `tumor_id` (identification number for the tumor)
* `record_time` (whether or not to record time points).

Within each simulation, there are 10,000 possible antigenicities, regularly spaced between 0 and 1. This model was used to generate the results shown in Figure 3.

Several outcomes are possible in this model, each of which is saved in particular directory.


1. Neither an adenoma or carcnoma evolves, in which cases results are saved in "no_tumors/"
2. A carcinoma evolved and exsisted for 1 year. These results are saved in "extant_carcinomas/"
3. A carcinoma evolved, but was removed by the immune. These results are saved in "extinct_carcinomas/"
4. A carcinoma evoloved with less than 1 year left in the simulation, and was not elimated. These results are saved in "controlled_carcinomas/"
5. An adenoma evolved, but was eliminated by the immune system before evoloving into a carcinoma. These results are saved in "extinct_adinomas/"
6. An adenoma evolved and existed for the remaineder of the simulation without evolving into a carcionm. These results are saved in "controlled_adenomas/"
7. Neither an adenoma or carcinoma formed, and the epithelial tissue became so antigenic that it was removed by immune predation. These results are saved in "extinct_epithelial/"

The results of each simultation are saved to csv files. A summary file is always created, and ends in "_summary.csv". Within this file, each row contains informaiton about each "species", where 0 = epithelial, 1=adenoma, 2=carcinoma. The first few columns contain the parameter values used in the simulation. Acutal results are recorded starting with the "species" column, which indicates the cell type:

* **origin**  indicates at whafounder_protectiont time that species arose
* **Time** is when the species got above a certain size detection threshold. The default size threshold is 100.
* **w_mean_ant** is the weighted mean antigenicity of the species, where the weights are the frequency of each phenotype in the species
* **w_mean_suppression** is the weighted mean immunesuppression of the species
* **w_mean_blockade** is the weighted mean immune blockade of the species;
* **mean_size** is the average size of each phenotype within the species at the end of the simulation
* **Total_Size** is the total size of eaach species at the end of the simulation
* **Total_Suppressive** is the total number of immunosuppressive cells associated with the species
* **Total_Blockade** is the total number of cells that can have the Blockade strategy
* **Total_Tcell** is the total number of T cells associated with each species
* **founder_ant** is the antigenicity of the phenotype that cell that founded the species, e.g. the first carcinoma cell;
* **founder_suppression** the amount of immune suppression of the founder cell
* **founder_protection** the amount of immune blockade of the founder cell.
* **n_pheno** the number of unique phenotypes present at the end of the simulation.

Note that **Total_Suppressive**, **Total_Blockade**, and **Total_Tcell** can be floats becuase they are based on `immunosuppression` and `blockade_protection`, which are between 0 and 1.

If the argument `record_time` is set to `True`, then similar information is recorded for each phenotype at a regular intervals. These results are saved in a file with the same name as the summary file, but without "_summary" in the name. Similar to the summary file, each row is a phenotype (as opposed to a species), and the first few columns contain information about the parameters used in the model. Phenotype specific values are recorded starting with the **blockade** column:

* **blockade** the strength of immune blockade. This is either 0 or `blockade_protection`, which is also recorded in the **blockade_param** column
* **suppression** the strength of immune suppression. This is either 0 or `immunosuppression`, which is also recorded in the **suppression_param** column
* **antigenicity**  how antigenic the phenotype is *before* accounting for immune suppression
* **n_drivers** number of driver mutations. Determines species
* **species** species of phenotype. 0=epithelial, 1=adenoma, 2=carcinoma
* **antigen_idx** an index to look up the phenotype's antigeniciy
* **parent** the id of phenotype that created this one.
* **t_cell_infiltration** strength of immune kill *after* accounting for immune suppression
* **div_rate** division rate *before* accounting for the incrased growth from recruiting immunosuppressive cells
* **net_r** net growth rate, i.e. division rate - (death + immune kill) *after* accounting for immune suppression
* **origin** the time at which the phenotype was created

Simulation times vary depending on parameter combinations, completeing within 3-20 minutes.

Example run:

    /anaconda3/bin/python3.7 mutation_model.py -dst_dir mutation_model_results/0.7p_0.7s -f_prefix 0.7p_0.7s -bp_mutation_rate 2.91e-09 -blockade_protection 0.7 -immunosuppression 0.7 -tumor_id 0 -record_time True

Results are saved to "mutation_model_results/0.7p_0.7s/"

The file "make_mutation_model_args.py" generates parameter combinations and writes each to a line in a bash file that can be used to run each simulation serially, or as an array job on a cluster. The variable that determines the size of the parameter grid is `N_VALS_FOR_EACH_STRATEGY`.
