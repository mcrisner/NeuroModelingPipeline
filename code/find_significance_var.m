function find_significance_var(inputfile, outputfile, lme)

% Set default value to not doing LME models
if nargin < 3
    lme = false;
end

% Load model data
models = load(inputfile).models;
main_stats = models.main_stats;
all_stats = models.all_stats;

% Row names for Excel output
colNames = {'Outcome', 'Predictor', 'β (std)', ...
            'β CI (std)', 'β p-value (raw)', 'β p-value (adj)', ...
            'Model p-value', 'Intercept p-value',  'β (raw)', 'β CI (raw)'};

varTypes = {'string', 'string', 'double', 'double', ...
            'double', 'double', 'double', 'string', 'double', 'string'};

% Get the field names (response variables) of the struct
resp_vars = fieldnames(main_stats);
pred_vars = fieldnames(main_stats.(resp_vars{1}));

% Define covariate variables
covariate_tbl = all_stats.(resp_vars{1}).(pred_vars{1});
covariates_of_interest = covariate_tbl.Properties.RowNames;
binary_vars = covariates_of_interest(isnan(covariate_tbl.std_beta));

significant_result_count = 0;

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

% Iterate over each outcome
for resp_idx = 1:numel(resp_vars)

    % Iterate over each predictor variable
    for var_idx = 1:length(pred_vars)
        mdl_table = main_stats.(resp_vars{resp_idx}).(pred_vars{var_idx});
    
        % Check if beta, model, and intercept p-values are significant (p <= 0.05)
        if mdl_table.model_pvalue(1) <= 0.05 && mdl_table.pvalue(1) <= 0.05 && mdl_table.Intercept_pvalue(1) <= 0.05
            significant_result_count = significant_result_count + 1;
    
            % Create string of confidence intervals
            std_beta_CI = ['[', num2str(round(mdl_table.std_beta_CI_low(1), 2)), ...
                           ', ', num2str(round(mdl_table.std_beta_CI_up(1), 2)), ']'];
            beta_CI = ['[', num2str(round(mdl_table.beta_CI_low(1), 2)), ', ', ...
                       num2str(round(mdl_table.beta_CI_up(1), 2)), ']'];
    
            % Create new spreadsheet row
            new_row = {mdl_table.Properties.RowNames{1}, pred_vars{var_idx}, ...
                       mdl_table.std_beta(1), std_beta_CI, ...
                       mdl_table.pvalue(1), mdl_table.qvalue(1), ...
                       mdl_table.model_pvalue(1), mdl_table.Intercept_pvalue(1), ...
                       mdl_table.regression_coef(1), beta_CI};
    
            adj_tbl = all_stats.(resp_vars{resp_idx}).(pred_vars{var_idx});
            adj_stats = cell(1, 4);
    
            for adj_idx = 1:numel(adj_tbl.Properties.RowNames)
    
                % Get adjustment beta if significant
                if adj_tbl.qvalue(adj_idx) <= 0.05
    
                    % Add name of covariate
                    adj_stats{1} = [adj_stats{1} adj_tbl.Properties.RowNames{adj_idx} ', '];
    
                    % Add adjusted p-value
                    qval = adj_tbl.qvalue(adj_idx);
                    if qval < 0.0005
                        adj_stats{2} = [adj_stats{2} num2str(sprintf('%.2e', qval)) ', '];
                    else
                        adj_stats{2} = [adj_stats{2} num2str(round(qval, 4)) ', '];
                    end
    
                    if ~ismember(adj_tbl.Properties.RowNames{adj_idx}, binary_vars)
                        % Add standardized beta
                        adj_stats{3} = [adj_stats{3} num2str(round(adj_tbl.std_beta(adj_idx), 4)) ', '];
    
                        % Add standardized CI
                        adj_stats{4} = [adj_stats{4} '[', num2str(round(adj_tbl.std_beta_CI_low(adj_idx), 4)), ...
                            ', ', num2str(round(adj_tbl.std_beta_CI_up(adj_idx), 4)), '], '];
                    else
                        % Add raw beta
                        adj_stats{3} = [adj_stats{3} num2str(round(adj_tbl.regression_coef(adj_idx), 4)) ', '];
    
                        % Add raw CI
                        adj_stats{4} = [adj_stats{4} '[', num2str(round(adj_tbl.beta_CI_low(adj_idx), 4)), ...
                            ', ', num2str(round(adj_tbl.beta_CI_up(adj_idx), 4)), '], '];
                    end
                end
            end
    
            % Remove trailing commas if present
            if ~isempty(adj_stats{1})
                adj_stats{1} = adj_stats{1}(1:end-2);
                adj_stats{2} = adj_stats{2}(1:end-2);
                adj_stats{3} = adj_stats{3}(1:end-2);
                adj_stats{4} = adj_stats{4}(1:end-2);
    
                % Add empty strings if cell are empty
            else
                adj_stats{1} = '';
                adj_stats{2} = '';
                adj_stats{3} = '';
                adj_stats{4} = '';
            end
    
            % Append covariate cells to the new spreadsheet row
            new_row = [new_row adj_stats];
    
            % Append mixed effect statistics if applicable
            if lme
                new_row = [new_row {mdl_table.lme_intercept_sig, ...
                           mdl_table.lme_intercept_CI, mdl_table.lme_slope_sig, ...
                           mdl_table.lme_slope_CI}];
            end
    
            % Add row to significant results table
            func_results = [func_results; new_row];
        end            
    end
end

% Save Excel file
output_dir = fileparts(outputfile);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
writetable(func_results, outputfile);
disp([num2str(significant_result_count), ' significant results.']);
disp(['Saved to ' outputfile]);

end