# NeuroModelingPipeline

This repository contains MATLAB code developed as part of a research study conducted at **Harvard Medical School** and **Boston Children’s Hospital**, funded by the **National Science Foundation**. The study investigates the neural correlates of social isolation in youth, leveraging morphometric features from structural MRI and topological features derived from resting-state fMRI using advanced methods from graph and network theory [(Rubinov & Sporns, 2010)](https://pubmed.ncbi.nlm.nih.gov/19819337/).

The modeling pipeline performs statistical analysis of brain-behavior associations at the connectome-, network-, node-, and structural-level. It supports flexible model configurations with both linear and mixed-effects modeling options, suitable for large-scale neuroimaging datasets such as the Adolescent Brain Cognitive Development (ABCD) Study.

This pipeline is designed to operate within a high-performance computing environment and utilizes **batch job submission to Harvard's O2 cluster via the SLURM job scheduler**. It includes options for customizing memory, walltime, and parallel execution using MATLAB's `parcluster` and `batch` functions.

All core modeling code is provided to support **open science** and facilitate **reproducible analyses** in the neuroscience community, with the goal of **accelerating discovery in brain-behavior research**.

A preprint of the study can be found [here](https://www.biorxiv.org/content/10.1101/2025.03.17.643680v1.full).

# `stat_modeling_main.m`

## Setup Instructions

To submit this function, either:

1. Configure and execute `run_stat_modeling_main.m`, or  
2. Instantiate a `parcluster` object with the following:

### 1. Memory Allocation
- Allocate at least 3 GB of memory per CPU.
- While 3 GB has yet to result in any crashes, consider using 4 GB for increased reliability.

### 2. Walltime
- Estimate the walltime based on the number of predictors of interest:
  - 1–2 hours per predictor of interest.

After the `parcluster` is initialized, define **ALL** input parameters and submit the job.

---

## Example Setup

```matlab
c = parcluster;
c.AdditionalProperties.MemPerCPU = '3000';
c.AdditionalProperties.WallTime = '04:00:00';
c.AdditionalProperties.QueueName = 'priority';

path_to_data = '...MATTHEW/modeling/modeling_code/data';
path_to_output = '...MATTHEW/modeling/modeling_code/results';
[define ALL other input parameters]

c.batch(@stat_modeling_main, 0, {path_to_data, path_to_output, ...});
```

## Input Parameters (Required)

> **NOTE**: The order of these parameters is very important. Please submit **all parameters** in the order listed below.

### `path_to_data` (string)
- Path to the folder where the necessary data files are stored.  
  **Example**:  
  `/n/data1/bch/radiology/stamoulis/COOPERATION_STUDY/NEW_RAs/MATTHEW/modeling/modeling_code/data`

### `path_to_output` (string)
- Path to the folder where the results will be saved.  
  **Example**:  
  `/n/data1/bch/radiology/stamoulis/COOPERATION_STUDY/NEW_RAs/MATTHEW/modeling/modeling_code/results`

### `main_dataset_file` (string)
- Filename for the main data table containing any independent variables.  
- This includes subject demographics, CBCL survey data, and scan parameters.  
- The file should contain a column of subject IDs named `src_subject_id`.  
  **Example**:  
  `scan_table_Y2.mat`

### `connectome_data_file` (string)
- Filename for the data table containing connectome-level brain properties.  
  **Example**:  
  `connectome_props_Y2.mat`

### `network_data_file` (string)
- Filename for the data table containing network-level brain properties.  
  **Example**:  
  `network_props_Y2.mat`

### `node_data_file` (string)
- Filename for the data table containing node-level brain properties.  
  **Example**:  
  `node_props_Y2.mat`

### `struc_data_file` (string)
- Filename for the data structure containing structural brain properties.  
  **Example**:  
  `ses_2Year_brain_struct_11-Apr-2023.mat`

### `roi_groups_file` (string)
- Filename for the data table containing network groupings of individual nodes.  
- This file associates each brain node with the resting-state network it belongs to.  
- The networks referenced are based on Yeo et al., 2011.  
  **Example**:  
  `ROIGroups_01-Aug-2024.mat`

### `sub_list_file` (string)
- Filename for the list of subject IDs to include in models.  
  **Example**:  
  `sub_list_year2_2809.mat`

### `basic_network_indices` (row vector of integers)
- Indices for basic networks to analyze from the `roi_groups_file`.  
- This vector specifies a subset of networks for nodal-level FDR correction and labeling.  
  **Example**:  
  `[1:9, 13:21, 26:32, 35:37]`

### `output_title` (string)
- Label used to name the output files.  
  **Example**:  
  `'social_withdrawal'` → generates files like `struc_results_social_withdrawal.xlsx`.

### `run_connectome` (0 or 1)
- Run models that predict connectome properties.  
  **1** to enable, **0** to disable.

### `run_network` (0 or 1)
- Run models that predict network properties.  
  **1** to enable, **0** to disable.

### `run_node` (0 or 1)
- Run models that predict node properties.  
  **1** to enable, **0** to disable.

### `run_structural` (0 or 1)
- Run models that predict structural properties.  
  **1** to enable, **0** to disable.

### `run_precovid_cohort` (0 or 1)
- Include only participants scanned before COVID.  
  **1** to enable, **0** to disable.

### `adjustment_vars` (cell array of strings)
- Additional adjustment variables to include in the model.
- Use the following **shorthand** codes:

### `adjustment_vars` (cell array of strings)

- Specifies additional variables to include in the model as adjustments.
- Each variable examined thus far has a shorthand abbreviation, automatically incorporated into both file names and spreadsheet tabs.
- These shorthand names should be used to define `adjustment_vars`:
	* 'preferalone' → would rather be alone than with others
	* 'int' → internalization t-score
	* 'ext' → externalization t-score
	* 'cohesion' → family cohesion
	* 'conflict' → family conflict
	* 'mh' → parental mental health
	* 'psb' → prosocial behavior
	* 'dim' → discrimination
	* 'anh' → anhedonia
	* 'fearsocial' → fear of social situations
	* 'esteem' → low self-esteem
	* 'embarrass' → self-conscious or easily embarrassed
	* 'bully' → problems with bullying
	* 'dep' → depression
	* 'anx' → anxiety
	* 'energized' → energized by large crowds
	* 'shy' → is shy
	* 'meet' → dislikes meeting new people
	* 'shymeet' → feels shy about meeting new people
	* friendgroup → has a regular friend group
- Set to an empty cell array `{}` to run the base model.  
  **Example**: `{'int', 'ext'}`

### `interactions` (2D cell array of strings)
- Specifies interactions between variables to include in the model.
- Use the shorthand names described under `adjustment_vars`.
- Write each interaction as a 1x2 cell array with a pair of variables.
- Interactions are included only if both variables are already present in the model.
  **Example**: `{{'preferalone', 'int'}, {'preferalone', 'ext'}}` or `{}`
  
### `base_controls` (cell array of strings)
- Base set of control variables to include in each model.
- `'percent_censored'` is automatically renamed depending on the selected run.
- Structural models automatically ignore `'percent_censored'` and `'scan_hour'`.
**Example**:  
`{'propensity', 'age', 'sex', 'race_ethnicity_twoCat', 'familyincome', 'bmi_zscore_sex', 'pds_score', 'percent_censored', 'scan_hour'}`

### `pred_of_interest` (cell array of strings)
- Predictor(s) of interest to be included in model permutations.
- Use the full variable name from `scan_table_Y2.mat`.  
  **Example**: `{'cbcl_prefer_alone', 'ksads_anhedonia'}`

### `scan_run` (string)
- Indicates which run to analyze: `'best'`, `'second'`, or `'both'`.
- If `'both'`, an "intersection" spreadsheet will be created.
  **Example**: `'both'`

### `lme` (logical)
- Enables or disables mixed-effects modeling (random slope + intercept).  
  **Example**: `true` or `false`

### `mixed_effect_var` (string)
- Variable used to group random effects.  
  **Example**: `'census_division_code'`

### `logistic` (logical)
- Whether to use logistic regression instead of linear modeling.  
  **Example**: `true` or `false`

### `thresholding_method` (string)
- Method used for thresholding:
  - `'median'`
  - `'percentile75'`
  - `'moderate_outlier'`
  - `'extreme'`
**Example**: `'moderate_outlier'`

### `num_cores` (integer)
- Number of CPU cores to allocate for parallelized node-level modeling.  
  **Example**: `4`

---

## `create_tables.m`

## Input Parameters (Required)

### `output_dir` (string)

- Path to the directory where the results will be saved.
- Ideally, this should match `path_to_data` in `stat_modeling_main`.


## Return Variables

### `scan_table`

- A table containing participant information from various surveys and data files.

### `connectome_props`

- A structure containing connectome-level properties for both best & second-best runs across each thresholding method.

### `network_props`

- A structure containing network-level properties for both best & second-best runs across each thresholding method.

### `node_props`

- A structure containing node-level properties for both best & second-best runs across each thresholding method.

> **Note**: Although these variables are returned, the primary purpose of this function is to **save** the output objects to the folder specified in `path_to_data` in `stat_modeling_main`.

> This file only needs to be run once, unless the table generation scripts are changed.
