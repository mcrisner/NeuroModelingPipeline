======================================================================
stat_modeling_main.m
======================================================================
SETUP INSTRUCTIONS:

To submit this function, either (1) configure and execute run_stat_modeling_main.m or 
				(2) instantiate a parcluster object with the following:

1. Memory Allocation
   - Allocate at least 3 GB of memory per CPU.
   - While 3 GB has yet to result in any crashes, consider using 4 GB for increased reliability.

2. Walltime
   - Estimate the walltime based on the number of predictors of interest:
      * 1-2 hours per predictor of interest.
   - Additional time may be required if the user's fair share value is low.

After the parcluster is initialized, define ALL input parameters and submit the job.

======================================================================
EXAMPLE SETUP:

c = parcluster;
c.AdditionalProperties.MemPerCPU = '3000';
c.AdditionalProperties.WallTime = '04:00:00';
c.AdditionalProperties.QueueName = 'priority';

path_to_data = '...MATTHEW/modeling/modeling_code/data';
path_to_output = '...MATTHEW/modeling/modeling_code/results';
[define ALL other input parameters]

c.batch(@stat_modeling_main, 0, {path_to_data, path_to_output, ...});

======================================================================
PARAMETERS (INPUT) - **MUST SUBMIT THESE WITH EACH JOB**

NOTE: The order of these parameters is very important. Please submit ALL parameters in the order that they appear here:

path_to_data (string)
   * Path to the folder where the necessary data files are stored.
   * Example: '/n/data1/bch/radiology/stamoulis/COOPERATION_STUDY/NEW_RAs/MATTHEW/modeling/modeling_code/data'

path_to_output (string)
   * Path to the folder where the results will be saved.
   * Example: '/n/data1/bch/radiology/stamoulis/COOPERATION_STUDY/NEW_RAs/MATTHEW/modeling/modeling_code/results'

main_dataset_file (string)
   * Filename for the main data table containing any independent variables.
   * This includes subject demographics, CBCL survey data, and scan parameters.
   * The file should contain a column of subject IDs named 'src_subject_id'.
   * Example: 'scan_table_Y2.mat'

connectome_data_file (string)
   * Filename for the data table containing connectome-level brain properties
   * Example: 'connectome_props_Y2.mat'

network_data_file (string)
   * Filename for the data table containing network-level brain properties
   * Example: 'network_props_Y2.mat'

node_data_file (string)
   * Filename for the data table containing node-level brain properties
   * Example: 'node_props_Y2.mat'

struc_data_file (string)
   * Filename for the data structure containing structural brain properties
   * Example: 'ses_2Year_brain_struct_11-Apr-2023.mat'

roi_groups_file (string)
   * Filename for the data table containing network groupings of individual nodes.
   * This file associates each brain node with the resting-state network it belongs to.
   * The networks referenced are based on Yeo et al., 2011.
   * Example: 'ROIGroups_01-Aug-2024.mat';

sub_list_file (string)
   * Filename for the list of subject IDs to include in models
   * Example: 'sub_list_year2_2809.mat'

basic_network_indices (row vector of integers)
   * Indices for basic networks to analyze from the roi_groups_file.
   * This vector specifies a subset of networks for nodal-level FDR correction and labeling.
   * Example: [1:9, 13:21, 26:32, 35:37]
   * The above example will extract rows 1-9, 13-21, etc. from roi_groups_file
   * for the sake of nodal FDR correction and labeling.

sub_list_file (string)
   * Filename for the list of subject IDs to include in models
   * Example: 'sub_list_year2_2809.mat'

output_title (string)
   * Label used to name the output files. For example, 'social_withdrawal' will generate output files like struc_results_social_withdrawal.xlsx.

run_connectome (0 or 1)
   * Indicates whether to run models that predict connectome properties.
   * 1 to enable, 0 to disable.

run_network (0 or 1)
   * Indicates whether to run models that predict network properties.
   * 1 to enable, 0 to disable.

run_node (0 or 1)
   * Indicates whether to run models that predict node properties.
   * 1 to enable, 0 to disable.

run_structural (0 or 1)
   * Indicates whether to run models that predict structural properties.
   * 1 to enable, 0 to disable.

run_precovid_cohort (0 or 1)
   * Indicates whether to include model permutations with only subjects scanned before COVID.
   * 1 to enable, 0 to disable.

adjustment_vars (cell array of strings)
   * Specifies additional variables to include in the model as adjustments.
   * Each variable examined thus far has a shorthand abbreviation, automatically incorporated into both file names and spreadsheet tabs.
   * These shorthand names should be used to define adjustment_vars:
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
   * Set to an empty cell array {} to run the base model.
   * Example: {'int', 'ext'} to run a model with the base controls, independent variable of interest, internalization, and externalization, or {} for just the base model.

interactions (2D cell array of strings)
   * Specifies interactions between variables to include in the model.
   * Use the shorthand names described under adjustment_vars.
   * Write each interaction as a 1x2 cell array with a pair of variables.
   * Interactions will be incorporated into a model if that model contains both variables, either as the predictor of interest or as adjustments.
   * Example: {{'preferalone', 'int'}, {'preferalone', 'ext'}} or {} for none.

base_controls (cell array of strings)
   * Base set of control variables to use with each model.
   * Including 'percent_censored' will automatically rename the variable based on whether a model is running best or second-best run.
   * Structural models will automatically ignore 'percent_censored' and 'scan_hour'.
   * Example: {'propensity', 'age', 'sex', 'race_ethnicity_twoCat', 'familyincome', 'bmi_zscore_sex', 'pds_score', 'percent_censored', 'scan_hour'}

pred_of_interest (cell array of strings)
   * Predictor(s) of interest to be included in the model permutations.
   * Use the full variable name as it appears in scan_table_Y2.mat
   * Example: {'cbcl_prefer_alone', 'ksads_anhedonia'}

scan_run (string)
   * Indicates which run to include in the analysis: 'best', 'second', or 'both'.
   * If set to 'both', an 'intersection' spreadsheet will be created containing associations significant in both runs.
   * Example: 'both'

lme (logical)
   * Whether to run models with random slope and intercept effects.
   * Example: true to enable, false to disable.

mixed_effect_var (string)
   * Variable used to group random effects.
   * Example: 'census_division_code'

logistic (logical)
   * Indicates whether to use logistic regression models.
   * Example: false to disable, true to enable.

thresholding_method (string)
   * Method used for thresholding: 'median', 'percentile75', 'moderate_outlier', or 'extreme'.
   * Example: 'moderate_outlier'

num_cores (integer)
   * Number of cores to allocate for parallelized nodal modeling jobs.
   * Example: 4

======================================================================
create_tables.m
======================================================================
PARAMETERS (INPUT) - **MUST SUBMIT THESE WITH EACH JOB**

output_dir (string)
   * Path to the directory where the results will be saved.
   * Ideally, this should be the same path as path_to_data in stat_modeling_main.

======================================================================
RETURN VARIABLES (OUTPUT)

scan_table
   * A table containing participant information from various surveys & data files.

connectome_props
   * A structure containing connectome-level properties for both best & second-best runs and each thresholding method.

network_props
   * A structure containing network-level properties for both best  & second-best runs and each thresholding method.

node_props
   * A structure containing node-level properties for both best & second-best runs and each thresholding method.

NOTE: Although the function also returns these variables, its primary role in the pipeline is to save these objects to the folder that path_to_data in stat_modeling_main points to.
NOTE: This file only needs to be run once, unless the table generation scripts are changed.
