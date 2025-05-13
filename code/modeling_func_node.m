function modeling_func_node(base_controls, pred_table, resp_table, sub_list, roi_groups, basic_network_indices, outputfile, vars_of_interest, adjustments, lme, mixed_effect_var, interaction, logistic)

% Description:
%   This function models the relationship between a set of predictors and various node-level 
%   properties, adding interaction terms, and including mixed effects. It generates model 
%   statistics for each property and saves them in a structured output file.
%
% Input Arguments:
%   base_controls     - Cell array of strings specifying base control variables.
%   pred_table        - Table containing predictor data with subject IDs and related variables.
%   resp_table        - Structure containing response (node-level properties) tables with subject IDs.
%   sub_list          - List of subject IDs to include in the models.
%   roi_groups        - Table defining ROI groups and their node indices for FDR correction.
%   outputfile        - String specifying the name of the output file (without extension).
%   pred_of_interest  - Cell array of strings specifying predictor(s) of interest.
%   adjustments       - Cell array of strings specifying additional covariates.
%   lme               - Logical flag to specify whether to include a mixed-effects model.
%   mixed_effect_var  - String specifying the variable for random effects if lme is true.
%   interaction       - Interaction cell array consisting of two variable names.
%   logistic          - Whether to run logistic regression instead of linear
%
% Process:
%   1. If specified, adds interaction terms between two variables.
%   2. Filters to relevant ROI groups for FDR correction.
%   3. Defines the list of node properties to model.
%   4. Combines puberty stages 4 and 5 into a single group.
%   5. Constructs model data by combining predictor and response tables.
%   6. For each predictor and node property, fits a linear or mixed-effects model.
%      - Calculates model statistics (e.g., regression coefficients, confidence intervals).
%      - Extracts p-values for FDR correction.
%   7. Applies FDR correction to p-values within specified ROI groups.
%   8. Saves the model statistics in a structured output file.
%
% Node Properties Modeled:
%   - degreeNodes
%   - clusterCoef_local
%   - centrality_eigenvector
%
% Dependencies:
%   Path to 'get_mdl_stats' function for extracting model statistics.
%
% Notes:
%   - If lme is true, the mixed effect variable is included in the model as a random effect.
%   - Saves the results to a directory called 'results/mat' in the current working directory.

% Add interaction terms to main table
if ~isempty(interaction)
    interaction_name = [interaction{1} '_x_' interaction{2}];
    adjustments = [adjustments interaction_name];
    pred_table.(interaction_name) = pred_table.(interaction{1}) .* pred_table.(interaction{2});
end

% Combine base controls and adjustments into a single array
controls = [base_controls adjustments];

% Filter to relevant ROI groups for FDR correction
basic_networks = roi_groups(basic_network_indices, :);

% Load table of brain properties & format subject IDs
prop_names = fieldnames(resp_table);
num_nodes = width(resp_table.(prop_names{1}));

% Make string of node names (e.g. node_1, node_2... node_1088)
node_vars = strings(num_nodes, 1);
for i = 1:num_nodes
    node_vars(i) = "node_" + i;
end

% Combine puberty stage 4 and 5 into one group (4)
pred_table.pds_score(pred_table.pds_score == 5) = 4;

% Extract predictor data from scan_table
pred_vars = [controls vars_of_interest];
if lme
    pred_vars = [pred_vars mixed_effect_var];
end
pred_data = pred_table(:, ['src_subject_id' pred_vars]);

% List of statistics to return
stat_lst = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
            'std_beta', 'std_beta_CI_low', 'std_beta_CI_up',...
            'SE', 'Wald', 'pvalue','Intercept_pvalue','model_pvalue'};
adj_stat_lst = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
                'std_beta', 'std_beta_CI_low', 'std_beta_CI_up',...
                'pvalue', 'model_pvalue'};

% Initialize structure to store results
prop_stats = struct;

% Initialize structure to store covariate statistics
all_stats = struct;

% Iterate over each predictor
for var_idx = 1:numel(vars_of_interest)
    
    % Iterate over each property
    for prop_idx = 1:numel(prop_names)

        disp(['Outcome: ' prop_names{prop_idx} '    Predictor: ' vars_of_interest{var_idx}]);

        % Extract node data from the prop structure
        prop_table = resp_table.(prop_names{prop_idx});
        prop_table.Properties.RowNames = strrep(prop_table.Properties.RowNames, 'sub-NDAR', 'NDAR_');
        prop_table.Properties.VariableNames = node_vars;
        prop_table.src_subject_id = prop_table.Properties.RowNames;

        % Initialize statistics table
        stat_tbl = array2table(NaN(num_nodes, length(stat_lst)), 'VariableNames', stat_lst, 'RowNames', node_vars);

        if lme
            stat_tbl.lme_intercept_sig = repmat("", num_nodes, 1);
            stat_tbl.lme_intercept_CI = repmat("", num_nodes, 1);
            stat_tbl.lme_slope_sig = repmat("", num_nodes, 1);
            stat_tbl.lme_slope_CI = repmat("", num_nodes, 1);
        end

        % Intialize temporary table for storing tables in parallel
        temp_tables = cell(num_nodes, 1);

        % Initialize list of p-values for FDR correction
        p_values = NaN(num_nodes, 1 + numel(controls));
    
        % Iterate over each node
        parfor node_idx = 1:numel(node_vars)

            % Intialize variables to avoid MATLAB warnings
            intercept_CI = NaN;
            slope_CI = NaN;
            intercept_significance = 'Significant';
            slope_significance = 'Significant';
            adj_stat_tbl = [];
            mdl = [];

            if mod(node_idx, 100) == 0
                disp(['    node: ' num2str(node_idx)]);
            end
            warning('OK');

            % Combine the predictor data and props data for creating the complete model table
            model_data = innerjoin(pred_data, prop_table(:, {'src_subject_id', node_vars{node_idx}}), 'Keys', 'src_subject_id');
            model_data = model_data(ismember(model_data.src_subject_id, sub_list), :);

            % Filter down to the provided list of subject IDs
            model_data.src_subject_id = [];

            % Construct the model formula
            mdl_spec = [node_vars{node_idx}, ' ~ ', strjoin(controls, ' + ') ' + ' vars_of_interest{var_idx}];
            
            if lme
                % Add random effects for intercept and slope
                mdl_spec = [mdl_spec ' + (' vars_of_interest{var_idx} ' | ' mixed_effect_var ')'];
    
                % Check for NaN values in response
                valid_indices = ~isnan(model_data.(node_vars{node_idx})) & ~isinf(model_data.(node_vars{node_idx}));
                model_data = model_data(valid_indices, :);

                % Fit the linear mixed effect model
                try
                    if logistic
                        mdl = fitglme(model_data, mdl_spec, 'Distribution', 'binomial', 'Link', 'logit');
                    else
                        mdl = fitlme(model_data, mdl_spec);
                    end
                catch me
                    warning(me.message)
                end
    
                % Get mixed effect statistics
                [~, ~, params] = covarianceParameters(mdl);
                intercept_lower = params{1}.Lower(1);
                intercept_upper = params{1}.Upper(1);
                intercept_CI = ['[' num2str(intercept_lower) ', ' num2str(intercept_upper) ']'];
                slope_lower = params{1}.Lower(3);
                slope_upper = params{1}.Upper(3);
                slope_CI = ['[' num2str(slope_lower) ', ' num2str(slope_upper) ']'];
    
                % Check significance of mixed effects
                if ((intercept_lower <= 0) && (intercept_upper >= 0)) || isnan(intercept_lower) || isnan(intercept_upper)
                    intercept_significance = 'Not significant';
                end
                if (slope_lower <= 0) && (slope_upper >= 0) || isnan(slope_lower) || isnan(slope_upper)
                    slope_significance = 'Not significant';
                end
            else
                % Fit the linear model
                try
                    if logistic
                        mdl = fitglm(model_data, mdl_spec, 'Distribution', 'binomial', 'Link', 'logit');
                    else
                        mdl = fitglm(model_data, mdl_spec);
                    end
                catch me
                    warning(me.message)
                end
            end
                 
            if strcmp(lastwarn, 'OK')
                % Initialize table to save adjustment statistics
                adj_stat_tbl = array2table(NaN(length(mdl.CoefficientNames)-2, length(adj_stat_lst)), 'VariableNames', adj_stat_lst, 'RowNames', mdl.CoefficientNames(2:end-1));
                stat_idx = find(strcmp(mdl.CoefficientNames, vars_of_interest{var_idx}));

                % Get statistics for predictor of interest and adjustments
                for coef = 2:length(mdl.CoefficientNames)   % Start at 2 to skip intercept
                    if coef == stat_idx
                        % Add model stats for current predictor to the main table
                        stat_row = get_mdl_stats(mdl, coef, stat_lst);
                        if lme
                            stat_row.lme_intercept_sig = intercept_significance;
                            stat_row.lme_intercept_CI = intercept_CI;
                            stat_row.lme_slope_sig = slope_significance;
                            stat_row.lme_slope_CI = slope_CI;
                        end
                        stat_tbl(node_idx, :) = stat_row;
                        p_values(node_idx, 1) = stat_row.pvalue;
                    else
                        % Add model stats for adjustment to the adjustments table
                        adj_stat_row = get_mdl_stats(mdl, coef, adj_stat_lst);
                        adj_stat_tbl(coef - 1, :) = adj_stat_row;
                    end
                end
            else
                % Save error message to log file
                log_file = 'model_errors.log';  % Define the log file name
                error_message = lastwarn;       % Get the last warning message
                timestamp = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');  % Current timestamp
                log_entry = sprintf('[%s] Predictor: %s    Outcome: %s    Warning: %s\n', ...
                                     timestamp, vars_of_interest{var_idx}, prop_names{prop_idx}, error_message);
            
                % Append the error message to the log file
                fid = fopen(log_file, 'a');  % Open the file in append mode
                if fid ~= -1
                    fprintf(fid, log_entry);
                    fclose(fid);
                else
                    warning('Unable to open the log file for writing.');
                end

                % Mock table to store NaNs for this model
                adj_stat_tbl = array2table(NaN(length(pred_data.Properties.VariableNames)-3, length(adj_stat_lst)), 'VariableNames', adj_stat_lst, 'RowNames', pred_data.Properties.VariableNames(2:end-2));
            end
            
            % Save adjustments table
            temp_tables{node_idx} = adj_stat_tbl;

        end

        % Add table to the main output structure
        prop_stats.(vars_of_interest{var_idx}).(prop_names{prop_idx}) = stat_tbl;

        % Assign models back to all_models structure outside the parfor loop
        for node_idx = 1:num_nodes
            tbl = temp_tables{node_idx};
            all_stats.(vars_of_interest{var_idx}).(prop_names{prop_idx}).(node_vars(node_idx)) = tbl;
            for row = 1:height(tbl)
                p_values(node_idx, row + 1) = tbl.pvalue(row);
            end
        end

        % Add a column of FDR corrected values
        for roi = 1:height(basic_networks)
            roi_nodes = basic_networks{roi, "Indices"};
            roi_p_values = p_values(cell2mat(roi_nodes), :);

            for var = 1:width(roi_p_values)
                q_values = mafdr(roi_p_values(:, var), 'BHFDR', true);

                for idx = 1:numel(roi_nodes{1})
                    node_idx = roi_nodes{1}(idx);
                    if var == 1
                        prop_stats.(vars_of_interest{var_idx}).(prop_names{prop_idx}).qvalue(node_idx) = q_values(idx);
                    else
                        all_stats.(vars_of_interest{var_idx}).(prop_names{prop_idx}).(node_vars(node_idx)).qvalue(var-1) = q_values(idx);
                    end
                end
            end
        end
    end
end

% Combine structures into a single output structure
models = struct;
models.prop_stats = prop_stats;
models.all_stats = all_stats;

% Save file
output_dir = fileparts(outputfile);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(outputfile, 'models');

end
