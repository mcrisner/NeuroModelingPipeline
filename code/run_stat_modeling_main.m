% ======================================================================================
% MODELING JOB SUBMISSION SCRIPT
%
% This script defines key parameters for the statistical modeling pipeline
% and submits a batch job for execution.

% ======================================================================================
% INITIALIZE PARALLEL COMPUTING CLUSTER

c = parcluster;                                  % Initialize parcluster object
c.AdditionalProperties.MemPerCPU = '3000';       % Allocate 3GB RAM per CPU
c.AdditionalProperties.WallTime  = '08:00:00';   % Set job time limit (8 hours)
c.AdditionalProperties.QueueName = 'priority';   % Assign job to priority queue

% ======================================================================================
% DEFINE INPUT PARAMETERS

% Path to the folder where the necessary data files are stored (string)
% Example: '/n/data1/bch/radiology/stamoulis/COOPERATION_STUDY/NEW_RAs/MATTHEW/modeling/modeling_code/data'
path_to_data          = '';

% Path to the folder where the results will be saved (string)
% Example: '/n/data1/bch/radiology/stamoulis/COOPERATION_STUDY/NEW_RAs/MATTHEW/modeling/modeling_code/results'
path_to_output        = '';

% Filename for the main data table containing any independent variables (string)
% e.g., subject demographics, CBCL survey data, scan parameters...
% Should include a column of subject IDs named 'src_subject_id'
% Example: 'scan_table_Y2.mat'
main_dataset_file     = '';

% Filenames for the data tables containing functional brain properties (strings)
% These files contain dependent variables related to different scales of functional brain organization:
% - Connectome-level (whole-brain connectivity)
% - Network-level (properties of large-scale networks)
% - Node-level (properties of individual nodes within networks)
% Example: 'connectome_props_Y2.mat'
connectome_data_file  = '';
network_data_file     = '';
node_data_file        = '';

% Filename for the data structure containing structural brain properties (string)
% This file contains dependent variables related to brain morphology.
% Example: 'ses_2Year_brain_struct_11-Apr-2023.mat'
struc_data_file       = '';

% Filename for the data table containing network membership of individual nodes (string)
% This file associates each brain node with the resting-state network it belongs to.
% The networks referenced are based on Yeo et al., 2011.
% Example: 'ROIGroups_01-Aug-2024.mat'
roi_groups_file       = '';

% Define indices for basic networks to analyze from the roi_groups_file (row vector of integers)
% This vector specifies a subset of networks for nodal-level FDR correction and labeling.
% Example: [1:9, 13:21, 26:32, 35:37]
basic_network_indices = [];

% Filename for the list of subject IDs to include in models (string)
% Example: 'sub_list_year2_2809.mat'
sub_list_file         = '';

% Output label (string - used for naming results files)
% Example: 'social_withdrawal' will generate output files like struc_results_social_withdrawal.xlsx.
output_title          = '';

% Flags to control which models to run (1 = enabled, 0 = disabled)
%
% The user has the option to run all analyses together (multi-scale functional connectome 
% and structural analyses) or select specific types of models.
% If only one analysis is desired, set the corresponding flag to 1 and the others to 0.
% If a model is not run, its associated data file can be set to an empty string ('').
run_connectome        = 0;        % Whether to run models that predict connectome properties
run_network           = 0;        % Whether to run models that predict network properties
run_node              = 0;        % Whether to run models that predict node properties
run_structural        = 1;        % Whether to run models that predict structural properties

% Base control variables (cell array of strings - included in every model)
% Percent censored and scan hour are automatically excluded from structural models.
% Including 'percent_censored' will automatically adapt to best vs. second best runs.
base_controls         = {};

% List of adjustment variables to include in all models (cell array of strings)
% Use empty cell array {} to run with only base_controls
% Recommended to use shorthand variable codes (e.g., {'int'} for 'cbcl_scr_syn_internal_t')
% The shorthand codes serve as suffixes for file naming purposes.
% See 'README' or 'VARIABLE MAPPING in stat_modeling_main' to find/edit shorthand codes
adjustment_vars       = {};

% List of variable interactions to include (2D cell array of strings - use empty {} to run base model)
% For interaction between var1 and var2: {{'var1', 'var2'}}
% If no interaction are tested, then set interacti
interactions          = {};       % See README or VARIABLE MAPPING in stat_modeling_main to find variable codes     

% Independent variables of interest (cell array of strings - independent variables)
% These are the primary predictors in the analysis. In this study, we examined two key 
% variables (being withdrawn and preferring solitude) but users can modify these to 
% any variable contained within main_dataset_file.mat.
vars_of_interest      = {};

% Additional modeling parameters
lme                   = true;                % Use linear mixed-effects models (LME)
mixed_effect_var      = '';                  % Variable for which random effects are included
logistic              = false;               % Use logistic regression? (False = linear regression)
thresholding_stat     = '';                  % Connectivity threshold: 'median', 'percentile75', 'moderate_outlier', or 'extreme'
num_cores             = 4;                   % Number of cores for parallelized node-level modeling
scan_run              = 'both';              % Scan selection: 'best', 'second', or 'both'
                                             % - If set to 'both', an 'intersection' spreadsheet will be 
                                             % - created containing associations significant in both runs.

% ======================================================================================
% SUBMIT MODELING JOB

% Job submission to Harvard's Parallel Computing Cluster using Slurm Workload Manager
j=c.batch(@stat_modeling_main, 0, {path_to_data, path_to_output, ...
    main_dataset_file, connectome_data_file, network_data_file, ...
    node_data_file, struc_data_file, roi_groups_file, basic_network_indices, ...
    output_title, run_connectome, run_network, run_node, run_structural, ...
    sub_list_file, adjustment_vars, interactions, ...
    base_controls, vars_of_interest, scan_run, lme, mixed_effect_var, ...
    logistic, thresholding_stat, num_cores});
    