function stat_modeling_main(path_to_data, path_to_output, main_dataset_file, ...
                            connectome_data_file, network_data_file, ...
                            node_data_file, struc_data_file, roi_groups_file, ...
                            basic_network_indices, output_title, ...
                            run_connectome, run_network, run_node, ...
                            run_structural, sub_list_file, adjustment_vars, interactions, ...
                            base_controls, vars_of_interest, scan_run, lme, ...
                            mixed_effect_var, logistic, thresholding_stat, ...
                            num_cores)
% ======================================================================================
% STAT MODELING MAIN
% 
% This script executes functional and structural statistical models to analyze brain properties.
% It includes the setup for model parameters, running batch jobs for functional and structural
% models, and processing results. Submit with run_stat_modeling_main.
%
% Main sections:
%   1. Parameter Definitions
%   2. Functional Modeling
%   3. Overlap Comparison
%   4. Structural Modeling

% ======================================================================================
% CREATE AND CONFIGURE CLUSTER

% Job submission to Harvard's Parallel Computing Cluster using Slurm Workload Manager
c=parcluster;
c.AdditionalProperties.QueueName='short';
c.AdditionalProperties.MemPerCPU='2000';
c.AdditionalProperties.WallTime='01:00:00';

% ======================================================================================
% VARIABLE MAPPING

% Define mapping between shorthand names and their corresponding variables.
% The shorthand names serve as suffixes for file naming purposes.
%
% INSTRUCTIONS FOR EDITING:
% - The 'shorthand' is a brief identifier used in file names.
% - The 'full variable name' is the corresponding full variable in pred_table.
% - To add a new variable, add its shorthand and corresponding full name to the end of each array.
variable_map = containers.Map( ...
    { ...   % SHORTHAND NAMES
    'preferalone', 'int', 'ext', 'cohesion', 'conflict', 'mh', 'psb', 'dim', ...
    'anh', 'fearsocial', 'esteem', 'embarrass', 'bully', 'dep', 'anx', ...
    'energized', 'shy', 'meet', 'shymeet', 'friendgroup'}, ...
    { ...   % FULL VARIABLE NAMES
        'cbcl_prefer_alone', ...
        'cbcl_scr_syn_internal_t', ...
        'cbcl_scr_syn_external_t', ...
        'fes_cohesion', ...
        'fes_conflict', ...
        'abcl_parent_mh', ...
        'psb_mean', ...
        'discrimination', ...
        'ksads_anhedonia', ...
        'ksads_fear_social', ...
        'ksads_low_esteem', ...
        'cbcl_self_conscious', ...
        'bullied', ...
        'cbcl_scr_dsm5_depress_t', ...
        'cbcl_scr_dsm5_anxdisord_t', ...
        'eatq_energized', ...
        'eatq_is_shy', ...
        'eatq_meet', ...
        'eatq_shy_meet', ...
        'ksads_friend_group'
    } ...
);
keys = variable_map.keys;

% Load control/covariate data
try
    pred_table = load(fullfile(path_to_data, main_dataset_file)).scan_table;
catch
    warning('Could not load control/covariate data.');
end

% Convert shorthand adjustment names to full variable names in pred_table
adjustments = {};
for var_idx = 1:numel(adjustment_vars)
    % If the provided adjustment is already a variable in pred_table,
    % add it to the array as is. Will take priority over variable map.
    if ismember(adjustment_vars{var_idx}, pred_table.Properties.VariableNames)
        adjustments = [adjustments adjustment_vars(var_idx)];
    else
        % Translate shorthand name based on variable_map
        for i = 1:length(keys)
            if strcmp(adjustment_vars{var_idx}, keys{i})
                adjustments = [adjustments {variable_map(keys{i})}];
            end
        end
    end
end

% Convert run_ variables to logicals
run_connectome = logical(run_connectome);
run_network = logical(run_network);
run_node = logical(run_node);
run_structural = logical(run_structural);

% Explicitely define list of runs to model
if strcmp(scan_run, 'best')
    run_list = {'best'};
elseif strcmp(scan_run, 'second')
    run_list = {'second'};
elseif strcmp(scan_run, 'both')
    run_list = {'best', 'second'};
end

% Replace strings in interactions with their corresponding values in variable_map
for i = 1:numel(interactions)
    for j = 1:numel(interactions{i})
        if isKey(variable_map, interactions{i}{j})
            interactions{i}{j} = variable_map(interactions{i}{j});
        else
            error('Key "%s" not found in variable_map.', interactions{i}{j});
        end
    end
end

% Load in list of subject IDs
try
    sub_list = load(fullfile(path_to_data, sub_list_file)).sub_list;
catch
    warning('Could not load list of subject IDs.');
end

% ======================================================================================
% RUN FUNCTIONAL MODELS

jobs = {};

% Load table of ROI groupings
try
    roi_groups = load(fullfile(path_to_data, roi_groups_file)).ROIGroups;
catch
    warning('Could not load ROI grouping data.');
end

for run_idx = 1:numel(run_list)

    run_name = run_list{run_idx};

    % Load connectome data
    try
        resp_table_connectome = load(fullfile(path_to_data, connectome_data_file)).connectome_props.(run_name).(thresholding_stat);
    catch
        warning('Could not load connectome data.');
    end

    % Load network data
    try
        resp_table_network = load(fullfile(path_to_data, network_data_file)).network_props.(run_name).(thresholding_stat);
    catch
        warning('Could not load network data.');
    end

    % Load node data
    try
        resp_table_node = load(fullfile(path_to_data, node_data_file)).node_props.(run_name).(thresholding_stat);
    catch
        warning('Could not load node data.');
    end
        
    % Modify 'percent_censored' control variable based on best/2ndbest run
    if strcmp(run_name, 'best')
        run_controls = strrep(base_controls, 'percent_censored', 'percent_censored_best');
    elseif strcmp(run_name, 'second')
        run_controls = strrep(base_controls, 'percent_censored', 'percent_censored_2ndbest');
    end

    % Create list of interactions to run with this covariate
    % This loop allows for running models with interactions betweenvariables
    % If no interactions are tested, loop will not be executed
    var_interactions = cell(1,1);
    for i = 1:numel(interactions)
        var1 = interactions{i}{1};
        var2 = interactions{i}{2};
        if ismember(var1, pred_table.Properties.VariableNames) && ismember(var2, pred_table.Properties.VariableNames)
            var_interactions = [var_interactions {interactions{i}}];
        end
    end

    outputfile_base = fullfile(path_to_output, 'mat', ['models_' output_title '_func_' run_name '_' strjoin(adjustment_vars, '_')]);    % Output file structure before being appended with e.g. '_connectome_lme.mat'

    % If no interactions are tested, loop will execute once
    for i = 1:numel(var_interactions)

        % Dynamically name output file
        if ~isempty(var_interactions{i})
            outputfile = [outputfile_base '_x'];
        else
            outputfile = outputfile_base;
        end

        % Run models
        if lme
            if run_connectome
                jobs{end+1} = c.batch(@modeling_func_connectome, 0, {run_controls, pred_table, resp_table_connectome, sub_list, [outputfile '_connectome_lme.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic});
            end
            if run_network
                jobs{end+1} = c.batch(@modeling_func_network, 0, {run_controls, pred_table, resp_table_network, sub_list, [outputfile '_network_lme.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic});
            end
            if run_node
                jobs{end+1} = c.batch(@modeling_func_node, 0, {run_controls, pred_table, resp_table_node, sub_list, roi_groups, basic_network_indices, [outputfile '_node_lme.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic}, 'Pool', num_cores);
            end
        else
            if run_connectome
                jobs{end+1} = c.batch(@modeling_func_connectome, 0, {run_controls, pred_table, resp_table_connectome, sub_list, [outputfile '_connectome.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic});
            end
            if run_network
                jobs{end+1} = c.batch(@modeling_func_network, 0, {run_controls, pred_table, resp_table_network, sub_list, [outputfile '_network.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic});
            end
            if run_node
                jobs{end+1} = c.batch(@modeling_func_node, 0, {run_controls, pred_table, resp_table_node, sub_list, roi_groups, basic_network_indices, [outputfile '_node.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic}, 'Pool', num_cores);
            end
        end
    end

    clearvars resp_table_connectome resp_table_network resp_table_node
    
end

for j = 1:length(jobs)
    wait(jobs{j});
end

% ======================================================================================
% FIND FUNCTIONAL RESULTS

jobs = {};

for run_idx = 1:numel(run_list)

    run_name = run_list{run_idx};

    % Create list of interactions to check for significance
    % This loop allows for running models with interactions betweenvariables
    % If no interactions are tested, loop will not be executed
    var_interactions = cell(1,1);
    for i = 1:numel(interactions)
        var1 = interactions{i}{1};
        var2 = interactions{i}{2};
        if ismember(var1, pred_table.Properties.VariableNames) && ismember(var2, pred_table.Properties.VariableNames)
            var_interactions = [var_interactions {interactions{i}}];
        end
    end

    inputfile_base = fullfile(path_to_output, 'mat', ['models_' output_title '_func_' run_name '_' strjoin(adjustment_vars, '_')]);
    outputfile = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_' run_name]);

    tabname_template = upper(strjoin(adjustment_vars, ' '));
    if isempty(tabname_template)
        tabname_template = 'BASE';
    end

    for i = 1:numel(var_interactions)

        % Dynamically name input file based on modeling output
        if ~isempty(var_interactions{i})
            inputfile = [inputfile_base '_x'];
        else
            inputfile = inputfile_base;
        end

        % Find significant results
        if run_connectome
            find_significance_func(inputfile, outputfile, roi_groups, basic_network_indices, 'connectome', lme, tabname_template);
        end
        if run_network
            find_significance_func(inputfile, outputfile, roi_groups, basic_network_indices, 'network', lme, tabname_template);
        end
        if run_node
            find_significance_func(inputfile, outputfile, roi_groups, basic_network_indices, 'node', lme, tabname_template)
        end
    end
end

for j = 1:length(jobs)
    wait(jobs{j});
end

% ======================================================================================
% OVERLAP COMPARISON

if strcmp(scan_run, 'both') || strcmp(scan_run, 'second')

    % Connectome-level
    if run_connectome
        best_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_best_connectome.xlsx']);
        second_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_second_connectome.xlsx']);
        output_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_intersection_connectome.xlsx']);
        compare_results(best_file, second_file, output_file, 'connectome');
    end

    % Network-level
    if run_network
        best_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_best_network.xlsx']);
        second_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_second_network.xlsx']);
        output_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_intersection_network.xlsx']);
        compare_results(best_file, second_file, output_file, 'network');
    end

    % Node-level
    if run_node
        best_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_best_node.xlsx']);
        second_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_second_node.xlsx']);
        output_file = fullfile(path_to_output, 'xlsx', ['func_results_' output_title '_intersection_node.xlsx']);
        compare_results(best_file, second_file, output_file, 'node');
    end

end

% ======================================================================================
% RUN STRUCTURAL MODELS

if run_structural

    struc_controls = setdiff(base_controls, {'percent_censored', 'scan_hour'}, 'stable');

    % Load structural data
    try
        struc_data = load(fullfile(path_to_data, struc_data_file));
    catch
        warning('Could not load structural data.');
    end

    jobs = {};

    % Create list of interactions to run with this covariate
    % This loop allows for running models with interactions betweenvariables
    % If no interactions are tested, loop will not be executed
    var_interactions = cell(1,1);
    for i = 1:numel(interactions)
        var1 = interactions{i}{1};
        var2 = interactions{i}{2};
        if ismember(var1, pred_table.Properties.VariableNames) && ismember(var2, pred_table.Properties.VariableNames)
            var_interactions = [var_interactions interactions(i)];
        end
    end

    % Name output file
    outputfile_base = fullfile(path_to_output, 'mat', ['models_' output_title '_struc_' strjoin(adjustment_vars, '_')]);

    for i = 1:numel(var_interactions)
        if ~isempty(var_interactions{i})
            outputfile = [outputfile_base '_x'];
        else
            outputfile = outputfile_base;
        end
        if lme
            jobs{end+1} = c.batch(@modeling_struc, 0, {struc_controls, pred_table, struc_data, sub_list, [outputfile '_lme.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic});
        else
            jobs{end+1} = c.batch(@modeling_struc, 0, {struc_controls, pred_table, struc_data, sub_list, [outputfile '.mat'], vars_of_interest, adjustments, lme, mixed_effect_var, var_interactions{i}, logistic});
        end
    end
    
    for j = 1:length(jobs)
        wait(jobs{j});
    end

end

% ======================================================================================
% FIND STRUCTURAL RESULTS

if run_structural
    
    % Create list of interactions to run with this covariate
    % This loop allows for running models with interactions betweenvariables
    % If no interactions are tested, loop will not be executed
    var_interactions = cell(1,1);
    for i = 1:numel(interactions)
        var1 = interactions{i}{1};
        var2 = interactions{i}{2};
        if ismember(var1, pred_table.Properties.VariableNames) && ismember(var2, pred_table.Properties.VariableNames)
            var_interactions = [var_interactions {interactions{i}}];
        end
    end

    inputfile_base = fullfile(path_to_output, 'mat', ['models_' output_title '_struc_' strjoin(adjustment_vars, '_')]);
    outputfile = fullfile(path_to_output, 'xlsx', ['struc_results_' output_title '.xlsx']);

    tabname_template = upper(strjoin(adjustment_vars, ' '));
    if isempty(tabname_template)
        tabname_template = 'BASE';
    end

    % Find significant results for each interaction term
    for i = 1:numel(var_interactions)
        if ~isempty(var_interactions{i})
            inputfile = [inputfile_base '_x'];
        else 
            inputfile = inputfile_base;
        end

        find_significance_struc(inputfile, outputfile, lme, tabname_template);
    end
end

end