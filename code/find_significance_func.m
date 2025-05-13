function find_significance_func(inputfile_template, outputfile_template, roi_groups, basic_network_indices, level_of_analysis, lme, tabname_template)

% Rename input/output files based on level of analysis (e.g. '_connectome.mat')
suffix = lower(level_of_analysis);
inputfile = sprintf('%s_%s%s', inputfile_template, suffix, '.mat');
outputfile = sprintf('%s_%s%s', outputfile_template, suffix, '.xlsx');

% Append '_lme' to input filename and modify tab name if using LME models
if lme
    inputfile = strrep(inputfile, '.mat', '_lme.mat');
    tabname = ['LME ' tabname_template];
else
    tabname = tabname_template;
end

% Load model data
models = load(inputfile).models;
prop_stats = models.prop_stats;
all_stats = models.all_stats;

% Row names for Excel output
colNames = {'Outcome', 'Predictor', 'Node', 'Network', 'β (std)', ...
            'β CI (std)', 'β p-value (raw)', 'β p-value (adj)', ...
            'Model p-value', 'Intercept p-value',  'β (raw)', 'β CI (raw)'};

varTypes = {'string', 'string', 'string', 'string', 'double', 'double', ...
            'double', 'double', 'double', 'string', 'double', 'string'};

% Get the field names (response variables) of the struct
pred_vars = fieldnames(prop_stats);
if strcmp(level_of_analysis, 'node')
    prop_names = fieldnames(prop_stats.(pred_vars{1}));
    node_names = prop_stats.(pred_vars{1}).(prop_names{1}).Properties.RowNames;
    covariate_tbl = all_stats.(pred_vars{1}).(prop_names{1}).(node_names{1});
elseif strcmp(level_of_analysis, 'network')
    network_names = fieldnames(prop_stats.(pred_vars{1}));
    prop_names = prop_stats.(pred_vars{1}).(network_names{1}).Properties.RowNames;
    covariate_tbl = all_stats.(pred_vars{1}).(network_names{1}).(prop_names{1});
    colNames(3) = [];
    varTypes(3) = [];
elseif strcmp(level_of_analysis, 'connectome')
    prop_names = prop_stats.(pred_vars{1}).Properties.RowNames;
    covariate_tbl = all_stats.(pred_vars{1}).(prop_names{1});
    colNames(3:4) = [];
    varTypes(3:4) = [];
end

% Define covariate and interaction variables
covariates_of_interest = covariate_tbl.Properties.RowNames;
binary_vars = covariates_of_interest(isnan(covariate_tbl.std_beta));
interaction_idx = find(contains(covariates_of_interest, '_x_'));

significant_result_count = 0;

% Add column names for interaction stats, if applicable
do_interaction = ~isempty(interaction_idx);
if do_interaction
    colNames = [colNames {'Interaction Term', 'Interaction β (std)', 'Interaction β CI (std)', ...
                          'Interaction p-value (adj)'}];
    varTypes = [varTypes {'string', 'string', 'string', 'string'}];
    tabname = [tabname ' X'];
end

% Add column names for covariates
colNames = [colNames {'Significant Covariates', 'Covariate p-values', ...
                      'Covariate βs (std)', 'Covariate β CIs (std)'}];
varTypes = [varTypes {'string', 'string', 'string', 'string'}];

% Add column names for mixed effect stats, if applicable
if lme
    colNames = [colNames {'Mixed Effect Intercept', 'Mixed Effect Intercept CI', ...
                          'Mixed Effect Slope', 'Mixed Effect Slope CI'}];
    varTypes = [varTypes {'string', 'string', 'string', 'string'}];
end

% Initialize table to store significant results
func_results = table('Size', [0 numel(varTypes)], 'VariableTypes', varTypes, 'VariableNames', colNames);

% Filter to relevant ROI groups
basic_networks = roi_groups(basic_network_indices, :);

% ==================== CONNECTOME ====================

if strcmp(level_of_analysis, 'connectome')

    % Iterate through each response variable
    for i = 1:length(pred_vars)
        mdl_table = prop_stats.(pred_vars{i});

        % Iterate through each predictor
        for prop_idx = 1:numel(prop_names)

            % Check if beta, model, and intercept p-values are significant (p <= 0.05)
            if ~strcmp(prop_names{prop_idx}, 'transitivity') && mdl_table.model_pvalue(prop_idx) <= 0.05 && mdl_table.qvalue(prop_idx) <= 0.05 && mdl_table.Intercept_pvalue(prop_idx) <= 0.05
                significant_result_count = significant_result_count + 1;

                % Create string of confidence intervals
                std_beta_CI = ['[', num2str(round(mdl_table.std_beta_CI_low(prop_idx), 4)), ...
                               ', ', num2str(round(mdl_table.std_beta_CI_up(prop_idx), 4)), ']'];
                beta_CI = ['[', num2str(round(mdl_table.beta_CI_low(prop_idx), 4)), ', ', ...
                           num2str(round(mdl_table.beta_CI_up(prop_idx), 4)), ']'];

                % Create new spreadsheet row
                new_row = {mdl_table.Properties.RowNames{prop_idx}, pred_vars{i}, ...
                           mdl_table.std_beta(prop_idx), std_beta_CI, ...
                           mdl_table.pvalue(prop_idx), mdl_table.qvalue(prop_idx), ...
                           mdl_table.model_pvalue(prop_idx), mdl_table.Intercept_pvalue(prop_idx), ...
                           mdl_table.regression_coef(prop_idx), beta_CI};

                % Extract table of covariates for this model
                covariate_tbl = all_stats.(pred_vars{i}).(mdl_table.Properties.RowNames{prop_idx});

                % Find the indices where row names match the target names
                [~, idx] = ismember(covariates_of_interest, covariate_tbl.Properties.RowNames);

                % Create cells that will later be appended to the new spreadsheet row
                covariate_stats = cell(1, 4);
                
                % Iterate over the covariates of interest
                for cov = 1:length(idx)
                    cov_idx = idx(cov);

                    % Extract whether a covariate was significant
                    if covariate_tbl.qvalue(cov_idx) <= 0.05
                        % Add name of covariate
                        covariate_stats{1} = [covariate_stats{1} covariates_of_interest{cov} ', '];
                        
                        % Add adjusted p-value
                        qval = covariate_tbl.qvalue(cov_idx);
                        if qval < 0.0005
                            covariate_stats{2} = [covariate_stats{2} num2str(sprintf('%.2e', qval)) ', '];
                        else
                            covariate_stats{2} = [covariate_stats{2} num2str(round(qval, 4)) ', '];
                        end
                        

                        if ~ismember(covariates_of_interest{cov}, binary_vars)
                            % Add standardized beta
                            covariate_stats{3} = [covariate_stats{3} num2str(round(covariate_tbl.std_beta(cov_idx), 4)) ', '];
    
                            % Add standardized CI
                            covariate_stats{4} = [covariate_stats{4} '[', num2str(round(covariate_tbl.std_beta_CI_low(cov_idx), 4)), ...
                                   ', ', num2str(round(covariate_tbl.std_beta_CI_up(cov_idx), 4)), '], '];
                        else
                            % Add raw beta
                            covariate_stats{3} = [covariate_stats{3} num2str(round(covariate_tbl.regression_coef(cov_idx), 4)) ', '];
    
                            % Add raw CI
                            covariate_stats{4} = [covariate_stats{4} '[', num2str(round(covariate_tbl.beta_CI_low(cov_idx), 4)), ...
                                   ', ', num2str(round(covariate_tbl.beta_CI_up(cov_idx), 4)), '], '];
                        end
                    end
                end
                
                % Remove trailing commas if present
                if ~isempty(covariate_stats{1})
                    covariate_stats{1} = covariate_stats{1}(1:end-2);
                    covariate_stats{2} = covariate_stats{2}(1:end-2);
                    covariate_stats{3} = covariate_stats{3}(1:end-2);
                    covariate_stats{4} = covariate_stats{4}(1:end-2);
                % Add empty strings if cell are empty
                else
                    covariate_stats{1} = '';
                    covariate_stats{2} = '';
                    covariate_stats{3} = '';
                    covariate_stats{4} = '';
                end

                % Gather statistics for the interaction term
                if do_interaction
                    interaction_stats = cell(1, 4);
                    interaction_stats{1} = 'Not Significant';
                    if covariate_tbl.qvalue(interaction_idx) <= 0.05
                        interaction_stats{1} = 'Significant';
                    end
                    interaction_stats{2} = num2str(round(covariate_tbl.std_beta(interaction_idx), 4));
                    interaction_stats{3} = ['[', num2str(round(covariate_tbl.std_beta_CI_low(interaction_idx), 4)), ...
                                            ', ', num2str(round(covariate_tbl.std_beta_CI_up(interaction_idx), 4)), ']'];
                    interaction_stats{4} = num2str(round(covariate_tbl.qvalue(interaction_idx), 4));

                    new_row = [new_row interaction_stats];
                end

                % Append covariate cells to the new spreadsheet row
                new_row = [new_row covariate_stats];

                % Append mixed effect statistics if applicable
                if lme
                    new_row = [new_row {mdl_table.lme_intercept_sig(prop_idx), ...
                               mdl_table.lme_intercept_CI(prop_idx), mdl_table.lme_slope_sig(prop_idx), ...
                               mdl_table.lme_slope_CI(prop_idx)}];
                end

                % Add row to significant results table
                func_results = [func_results; new_row];
            end            
        end
    end

% ==================== NETWORK ====================

elseif strcmp(level_of_analysis, 'network')

    % Iterate through each response variable
    for i = 1:length(pred_vars)
        pred_struct = prop_stats.(pred_vars{i});
        
        % Iterate through each predictor
        for prop_idx = 1:numel(prop_names)
    
            % Iterate through each network
            for n = 1:numel(network_names)
            
                mdl_table = pred_struct.(network_names{n});
    
                % Check if beta, model, and intercept p-values are significant (p <= 0.05)
                if ~strcmp(prop_names{prop_idx}, 'transitivity') && mdl_table.model_pvalue(prop_idx) <= 0.05 && mdl_table.qvalue(prop_idx) <= 0.05 && mdl_table.Intercept_pvalue(prop_idx) <= 0.05
                    significant_result_count = significant_result_count + 1;

                    % Create string of confidence intervals
                    std_beta_CI = ['[', num2str(round(mdl_table.std_beta_CI_low(prop_idx), 4)), ...
                                   ', ', num2str(round(mdl_table.std_beta_CI_up(prop_idx), 4)), ']'];
                    beta_CI = ['[', num2str(round(mdl_table.beta_CI_low(prop_idx), 4)), ', ', ...
                               num2str(round(mdl_table.beta_CI_up(prop_idx), 4)), ']'];

                    % Create new spreadsheet row
                    new_row = {mdl_table.Properties.RowNames{prop_idx}, pred_vars{i}, network_names{n}, ...
                               mdl_table.std_beta(prop_idx), std_beta_CI, ...
                               mdl_table.pvalue(prop_idx), mdl_table.qvalue(prop_idx), ...
                               mdl_table.model_pvalue(prop_idx), mdl_table.Intercept_pvalue(prop_idx), ...
                               mdl_table.regression_coef(prop_idx), beta_CI};
    
                    % Extract table of covariates for this model
                    covariate_tbl = all_stats.(pred_vars{i}).(network_names{n}).(mdl_table.Properties.RowNames{prop_idx});
    
                    % Find the indices where row names match the target names
                    [~, idx] = ismember(covariates_of_interest, covariate_tbl.Properties.RowNames);
    
                    % Create cells that will later be appended to the new spreadsheet row
                    covariate_stats = cell(1, 4);
                    
                    % Iterate over the covariates of interest
                    for cov = 1:length(idx)
                        cov_idx = idx(cov);
    
                        % Extract whether a covariate was significant
                        if covariate_tbl.qvalue(cov_idx) <= 0.05
                            % Add name of covariate
                            covariate_stats{1} = [covariate_stats{1} covariates_of_interest{cov} ', '];

                            % Add adjusted p-value
                            qval = covariate_tbl.qvalue(cov_idx);
                            if qval < 0.0005
                                covariate_stats{2} = [covariate_stats{2} num2str(sprintf('%.2e', qval)) ', '];
                            else
                                covariate_stats{2} = [covariate_stats{2} num2str(round(qval, 4)) ', '];
                            end
    
                            if ~ismember(covariates_of_interest{cov}, binary_vars)
                                % Add standardized beta
                                covariate_stats{3} = [covariate_stats{3} num2str(round(covariate_tbl.std_beta(cov_idx), 4)) ', '];
        
                                % Add standardized CI
                                covariate_stats{4} = [covariate_stats{4} '[', num2str(round(covariate_tbl.std_beta_CI_low(cov_idx), 4)), ...
                                       ', ', num2str(round(covariate_tbl.std_beta_CI_up(cov_idx), 4)), '], '];
                            else
                                % Add raw beta
                                covariate_stats{3} = [covariate_stats{3} num2str(round(covariate_tbl.regression_coef(cov_idx), 4)) ', '];
        
                                % Add raw CI
                                covariate_stats{4} = [covariate_stats{4} '[', num2str(round(covariate_tbl.beta_CI_low(cov_idx), 4)), ...
                                       ', ', num2str(round(covariate_tbl.beta_CI_up(cov_idx), 4)), '], '];
                            end
                        end
                    end
                    
                    % Remove trailing commas if present
                    if ~isempty(covariate_stats{1})
                        covariate_stats{1} = covariate_stats{1}(1:end-2);
                        covariate_stats{2} = covariate_stats{2}(1:end-2);
                        covariate_stats{3} = covariate_stats{3}(1:end-2);
                        covariate_stats{4} = covariate_stats{4}(1:end-2);
                    % Add empty strings if cell are empty
                    else
                        covariate_stats{1} = '';
                        covariate_stats{2} = '';
                        covariate_stats{3} = '';
                        covariate_stats{4} = '';
                    end

                    % Gather statistics for the interaction term
                    if do_interaction
                        interaction_stats = cell(1, 4);
                        interaction_stats{1} = 'Not Significant';
                        if covariate_tbl.qvalue(interaction_idx) <= 0.05
                            interaction_stats{1} = 'Significant';
                        end
                        interaction_stats{2} = num2str(round(covariate_tbl.std_beta(interaction_idx), 4));
                        interaction_stats{3} = ['[', num2str(round(covariate_tbl.std_beta_CI_low(interaction_idx), 4)), ...
                                                ', ', num2str(round(covariate_tbl.std_beta_CI_up(interaction_idx), 4)), ']'];
                        interaction_stats{4} = num2str(round(covariate_tbl.qvalue(interaction_idx), 4));
    
                        new_row = [new_row interaction_stats];
                    end
    
                    % Append covariate cells to the new spreadsheet row
                    new_row = [new_row covariate_stats];

                    % Append mixed effect statistics if applicable
                    if lme
                        new_row = [new_row {mdl_table.lme_intercept_sig(prop_idx), ...
                                   mdl_table.lme_intercept_CI(prop_idx), mdl_table.lme_slope_sig(prop_idx), ...
                                   mdl_table.lme_slope_CI(prop_idx)}];
                    end
    
                    % Add row to significant results table
                    func_results = [func_results; new_row];
                end
            end
        end
    end

% ==================== NODE ====================

elseif strcmp(level_of_analysis, 'node')
    
    % Iterate through each response variable
    for i = 1:length(pred_vars)
        pred_struct = prop_stats.(pred_vars{i});
        
        % Iterate through each predictor
        for prop_idx = 1:numel(prop_names)

            mdl_table = pred_struct.(prop_names{prop_idx});
    
            % Iterate through each network
            for n = 1:numel(node_names)
            
                % Check if beta, model, and intercept p-values are significant (p <= 0.05)
                if mdl_table.model_pvalue(n) <= 0.05 && mdl_table.qvalue(n) <= 0.05 && mdl_table.Intercept_pvalue(n) <= 0.05
                    significant_result_count = significant_result_count + 1;

                    % Create string of confidence intervals
                    std_beta_CI = ['[', num2str(round(mdl_table.std_beta_CI_low(n), 4)), ...
                                   ', ', num2str(round(mdl_table.std_beta_CI_up(n), 4)), ']'];
                    beta_CI = ['[', num2str(round(mdl_table.beta_CI_low(n), 4)), ', ', ...
                               num2str(round(mdl_table.beta_CI_up(n), 4)), ']'];

                    % Find network that node is part of
                    node_network = basic_networks.ROIName(cellfun(@(x) ismember(n, x), basic_networks.Indices));

                    % Create new spreadsheet row
                    new_row = {prop_names{prop_idx}, pred_vars{i}, num2str(n), node_network...
                               mdl_table.std_beta(n), std_beta_CI, ...
                               mdl_table.pvalue(n), mdl_table.qvalue(n), ...
                               mdl_table.model_pvalue(n), mdl_table.Intercept_pvalue(n), ...
                               mdl_table.regression_coef(n), beta_CI};
    
                    % Extract table of covariates for this model
                    covariate_tbl = all_stats.(pred_vars{i}).(prop_names{prop_idx}).(mdl_table.Properties.RowNames{n});
    
                    % Find the indices where row names match the target names
                    [~, idx] = ismember(covariates_of_interest, covariate_tbl.Properties.RowNames);
    
                    % Create cells that will later be appended to the new spreadsheet row
                    covariate_stats = cell(1, 4);
                    
                    % Iterate over the covariates of interest
                    for cov = 1:length(idx)
                        cov_idx = idx(cov);
    
                        % Extract whether a covariate was significant
                        if covariate_tbl.qvalue(cov_idx) <= 0.05
                            % Add name of covariate
                            covariate_stats{1} = [covariate_stats{1} covariates_of_interest{cov} ', '];

                            % Add adjusted p-value
                            qval = covariate_tbl.qvalue(cov_idx);
                            if qval < 0.0005
                                covariate_stats{2} = [covariate_stats{2} num2str(sprintf('%.2e', qval)) ', '];
                            else
                                covariate_stats{2} = [covariate_stats{2} num2str(round(qval, 4)) ', '];
                            end
    
                            if ~ismember(covariates_of_interest{cov}, binary_vars)
                                % Add standardized beta
                                covariate_stats{3} = [covariate_stats{3} num2str(round(covariate_tbl.std_beta(cov_idx), 4)) ', '];
        
                                % Add standardized CI
                                covariate_stats{4} = [covariate_stats{4} '[', num2str(round(covariate_tbl.std_beta_CI_low(cov_idx), 4)), ...
                                       ', ', num2str(round(covariate_tbl.std_beta_CI_up(cov_idx), 4)), '], '];
                            else
                                % Add raw beta
                                covariate_stats{3} = [covariate_stats{3} num2str(round(covariate_tbl.regression_coef(cov_idx), 4)) ', '];
        
                                % Add raw CI
                                covariate_stats{4} = [covariate_stats{4} '[', num2str(round(covariate_tbl.beta_CI_low(cov_idx), 4)), ...
                                       ', ', num2str(round(covariate_tbl.beta_CI_up(cov_idx), 4)), '], '];
                            end
                        end
                    end
                    
                    % Remove trailing commas if present
                    if ~isempty(covariate_stats{1})
                        covariate_stats{1} = covariate_stats{1}(1:end-2);
                        covariate_stats{2} = covariate_stats{2}(1:end-2);
                        covariate_stats{3} = covariate_stats{3}(1:end-2);
                        covariate_stats{4} = covariate_stats{4}(1:end-2);
                    % Add empty strings if cell are empty
                    else
                        covariate_stats{1} = '';
                        covariate_stats{2} = '';
                        covariate_stats{3} = '';
                        covariate_stats{4} = '';
                    end

                    % Gather statistics for the interaction term
                    if do_interaction
                        interaction_stats = cell(1, 4);
                        interaction_stats{1} = 'Not Significant';
                        if covariate_tbl.qvalue(interaction_idx) <= 0.05
                            interaction_stats{1} = 'Significant';
                        end
                        interaction_stats{2} = num2str(round(covariate_tbl.std_beta(interaction_idx), 4));
                        interaction_stats{3} = ['[', num2str(round(covariate_tbl.std_beta_CI_low(interaction_idx), 4)), ...
                                                ', ', num2str(round(covariate_tbl.std_beta_CI_up(interaction_idx), 4)), ']'];
                        interaction_stats{4} = num2str(round(covariate_tbl.qvalue(interaction_idx), 4));
    
                        new_row = [new_row interaction_stats];
                    end
    
                    % Append covariate cells to the new spreadsheet row
                    new_row = [new_row covariate_stats];

                    % Append mixed effect statistics if applicable
                    if lme
                        new_row = [new_row {mdl_table.lme_intercept_sig(n), ...
                                   mdl_table.lme_intercept_CI(n), mdl_table.lme_slope_sig(n), ...
                                   mdl_table.lme_slope_CI(n)}];
                    end
    
                    % Add row to significant results table
                    func_results = [func_results; new_row];
                end
            end
        end
    end
end

% Save Excel file
output_dir = fileparts(outputfile);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
writetable(func_results, outputfile, 'Sheet', tabname);
disp([num2str(significant_result_count), ' significant results.']);
disp(['Saved to ' outputfile]);

end