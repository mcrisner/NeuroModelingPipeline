function modeling_var_associations(base_controls, scan_table, outputfile, pred_of_interest, resp_vars, before_covid_only, lme, mixed_effect_var, logistic, combine_pds_score)

% Filter down to only subjects scanned before COVID onset (3/11/20), if applicable
if before_covid_only
    scan_table = scan_table(scan_table.after_pandemic == 0, :);
end

% Combine puberty stage 4 and 5 into one group (4)
if combine_pds_score
    scan_table.pds_score(scan_table.pds_score == 5) = 4;
end

% Extract predictor data from scan_table
pred_vars = [base_controls pred_of_interest];
if lme
    pred_vars = [pred_vars mixed_effect_var];
end
model_data = scan_table(:, [pred_vars resp_vars]);

% List of statistics to return
stat_lst = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
            'std_beta', 'std_beta_CI_low', 'std_beta_CI_up',...
            'pvalue','Intercept_pvalue','model_pvalue'};
adj_stat_lst = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
                'std_beta', 'std_beta_CI_low', 'std_beta_CI_up',...
                'pvalue','Intercept_pvalue','model_pvalue'};

% Initialize structure to store results
main_stats = struct;

% Initialize structure to store covariate statistics
all_stats = struct;

% Iterate over each outcome
for resp_idx = 1:numel(resp_vars)

    resp_var = resp_vars{resp_idx};

    p_values = NaN(numel(pred_of_interest), 1 + numel(base_controls));

    % Iterate over each predictor
    for var_idx = 1:numel(pred_of_interest)
    
        % Initialize statistics table
        stat_tbl = array2table(NaN(1, length(stat_lst)), 'VariableNames', stat_lst, 'RowNames', {resp_var});
       
        if lme
            stat_tbl.lme_intercept_sig = repmat("", height(stat_tbl), 1);
            stat_tbl.lme_intercept_CI = repmat("", height(stat_tbl), 1);
            stat_tbl.lme_slope_sig = repmat("", height(stat_tbl), 1);
            stat_tbl.lme_slope_CI = repmat("", height(stat_tbl), 1);
        end

        if contains(pred_of_interest{var_idx}, base_controls)
            var_controls = base_controls(~contains(base_controls,pred_of_interest{var_idx}));
        else
            var_controls = base_controls;
        end
        
        warning('OK');
    
        mdl_spec = [resp_var, ' ~ ', var_controls{1}];
        for i = 2:(numel(var_controls))    
            mdl_spec = [mdl_spec, ' + ', var_controls{i}];
        end
        mdl_spec = [mdl_spec, ' + ', pred_of_interest{var_idx}];
    
        if lme
            % Add random effect for intercept
            mdl_spec = [mdl_spec ' + (1 | ' mixed_effect_var ')'];
    
            % Add random effect for slope
            mdl_spec = [mdl_spec ' + (' pred_of_interest{var_idx} ' - 1 | ' mixed_effect_var ')'];
    
            % Check for NaN values in response
            valid_indices = ~isnan(model_data.(resp_var)) & ~isinf(model_data.(resp_var));
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
            intercept_lower = params{1}.Lower;
            intercept_upper = params{1}.Upper;
            intercept_CI = ['[' num2str(intercept_lower) ', ' num2str(intercept_upper) ']'];
            slope_lower = params{2}.Lower;
            slope_upper = params{2}.Upper;
            slope_CI = ['[' num2str(slope_lower) ', ' num2str(slope_upper) ']'];
    
            % Check significance of mixed effects
            intercept_significance = 'Significant';
            slope_significance = 'Significant';
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
            adj_stat_tbl = array2table(NaN(length(mdl.CoefficientNames)-2, length(adj_stat_lst)), 'VariableNames', adj_stat_lst, 'RowNames', var_controls);
            
            stat_idx = find(strcmp(mdl.CoefficientNames, pred_of_interest{var_idx}));
    
            % Get statistics for predictor of interest and adjustments
            for coef = 2:length(mdl.CoefficientNames)   % Start at 2 to skip intercept
                if coef == stat_idx
                    % Add model stats for current predictor to the main table
                    stat_tbl = get_mdl_stats(mdl, coef, stat_lst);
                    if lme
                        stat_tbl.lme_intercept_sig = intercept_significance;
                        stat_tbl.lme_intercept_CI = intercept_CI;
                        stat_tbl.lme_slope_sig = slope_significance;
                        stat_tbl.lme_slope_CI = slope_CI;
                    end
                    p_values(var_idx, 1) = stat_tbl.pvalue;
                else
                    % Add model stats for adjustment to the adjustments table
                    adj_stat_row = get_mdl_stats(mdl, coef, adj_stat_lst);
                    adj_stat_tbl(coef - 1, :) = adj_stat_row;
                    p_values(var_idx, coef) = adj_stat_row.pvalue;
                end
            end    
        else
            % Save error message to log file
            log_file = 'model_errors.log';  % Define the log file name
            error_message = lastwarn;       % Get the last warning message
            timestamp = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');  % Current timestamp
            log_entry = sprintf('[%s] Predictor: %s    Outcome: %s    Warning: %s\n', ...
                                 timestamp, pred_of_interest{var_idx}, resp_vars{resp_idx}, error_message);
            disp(error_message)
        
            % Append the error message to the log file
            fid = fopen(log_file, 'a');  % Open the file in append mode
            if fid ~= -1
                fprintf(fid, log_entry);
                fclose(fid);
            else
                warning('Unable to open the log file for writing.');
            end
    
            % Mock table to store NaNs for this model
            adj_stat_tbl = array2table(NaN(length(var_controls), length(adj_stat_lst)), 'VariableNames', adj_stat_lst, 'RowNames', var_controls);
        end
    
        % Add table to the main output structure
        all_stats.(resp_var).(pred_of_interest{var_idx}) = adj_stat_tbl;
    
        % Add table to the main output structure
        main_stats.(resp_var).(pred_of_interest{var_idx}) = stat_tbl;
    
    end

    % Add a column of FDR corrected values
    for var = 1:width(p_values)
        q_values = mafdr(p_values(:, var), 'BHFDR', true);
        
        for var_idx = 1:numel(pred_of_interest)
            if ~isnan(q_values(var_idx))
                if var == 1
                    main_stats.(resp_var).(pred_of_interest{var_idx}).qvalue = q_values(var_idx);
                else
                    all_stats.(resp_var).(pred_of_interest{var_idx}).qvalue(var-1) = q_values(var_idx);
                end
            end
        end
    end
end

% Combine structures into a single output structure
models = struct;
models.main_stats = main_stats;
models.all_stats = all_stats;

% Save file
output_dir = fileparts(outputfile);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(outputfile, 'models');

end