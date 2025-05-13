% =========================================================================
% VARIABLE ASSOCIATION MODELING SCRIPT
%
% This script loads the main data table and defines key variables for 
% statistical modeling of associations between predictor variables and 
% response variables. It then runs the modeling function and extracts 
% significant results for output to an Excel spreadsheet.

% =========================================================================
% DEFINE INPUT PARAMETERS

% Path of where to find input datasets
path_to_data = '';

% Load the main data table containing subject-level variables
scan_table = load(fullfile(path_to_data, ''));

% List of independent variables of interest 
% (Includes demographic, social, and behavioral factors)
vars_of_interest = {''};

% Response variables to model (social withdrawal-related measures)
resp_vars = {''};

% Set of base control variables (potential confounds adjusted for in every model)
base_controls = {''};

% =========================================================================
% MODEL CONFIGURATION OPTIONS

before_covid_only = false;        % Restrict analysis to pre-pandemic data only
lme = false;                      % Use linear mixed-effects models (LME)
mixed_effect_var = 'site';       % Variable for grouping random effects
logistic = false;                % Use logistic regression? (False = linear regression)
combine_pds_score = true;       % Combine PDS score categories? (False = separate)

% =========================================================================
% OUTPUT FILE LOCATIONS

% Define output file paths for saving results
outputfile_mat = fullfile(''); % Save model output
outputfile_xlsx = fullfile(''); % Save significant results

% =========================================================================
% RUN MODELING AND SIGNIFICANCE ANALYSIS

% Step 1: Run the variable association modeling
modeling_var_associations(base_controls, scan_table, outputfile_mat, ...
                          vars_of_interest, resp_vars, before_covid_only, ...
                          lme, mixed_effect_var, logistic, combine_pds_score);

% Step 2: Extract significant results and save to Excel
find_significance_var(outputfile_mat, outputfile_xlsx, lme);
