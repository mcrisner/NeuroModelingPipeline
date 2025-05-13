function find_significance_struc(inputfile_template, outputfile, lme, tabname_template)

% Append '_lme.mat' to input filename and modify tab name if using LME models
if lme
    inputfile = [inputfile_template '_lme.mat'];
    tabname = ['LME ' tabname_template];
else
    inputfile = [inputfile_template '.mat'];
    tabname = tabname_template;
end

% Load model data
models = load(inputfile).models;
prop_stats = models.prop_stats;
all_stats = models.all_stats;

% Get the field names (response variables) of the struct
pred_vars = fieldnames(prop_stats);
all_regions = fieldnames(prop_stats.(pred_vars{1}));
cort_regions = all_regions(1:68);
prop_names = prop_stats.(pred_vars{1}).(all_regions{1}).Properties.RowNames;

% Define covariate and interaction variables
covariate_tbl = all_stats.(pred_vars{1}).(all_regions{1}).(prop_names{1});
covariates_of_interest = covariate_tbl.Properties.RowNames;
binary_vars = covariates_of_interest(isnan(covariate_tbl.std_beta));
interaction_idx = find(contains(covariates_of_interest, '_x_'));

% Row names for Excel output
colNames = {'Region', 'Predictor', 'THICK: β', 'THICK: β CI', 'THICK: β q-val', ...
            'THICK: Sig. Cov.', 'THICK: Cov. q-vals', 'THICK: Cov. β', ...
            'THICK: Cov. β CI', 'THICK: Ixn. β', 'THICK: Ixn. β CI', ...
            'THICK: Ixn. q-val', 'VOL: β', 'VOL: β CI', ...
            'VOL: β q-val', 'VOL: Sig. Cov.' 'VOL: Cov. β', 'VOL: Cov. β CI', ...
            'VOL: Cov. q-vals', 'VOL: Ixn. β', 'VOL: Ixn. β CI', ...
            'VOL: Ixn. q-val', 'WM: β', 'WM: β CI', 'WM: β q-val', ...
            'WM: Sig. Cov.' 'WM: Cov. β', 'WM: Cov. β CI', ...
            'WM: Cov. q-val', 'WM: Ixn. β', 'WM: Ixn. β CI', ...
            'WM: Ixn. q-val'};
varTypes = {'string', 'string', 'double', 'string', 'double', 'string', 'string', ...
            'string', 'string', 'double', 'string', 'double', 'cell', ...
            'string', 'cell', 'string', 'string', 'string', 'string', ...
            'cell', 'string', 'cell', 'double', 'string', 'double', ...
            'string', 'string', 'string', 'string', 'cell', 'string', ...
            'cell'};

do_interaction = ~isempty(interaction_idx);
if do_interaction
    tabname = [tabname ' X'];
else
    ixn_indices = contains(colNames, 'Ixn.');
    colNames(ixn_indices) = [];
    varTypes(ixn_indices) = [];
end

% Dictionary for translating DK regions
decoded_regions = {
    'lh-banks-of-superior-temporal-sulcus', 'lh-caudalanteriorcingulate', ...
    'lh-caudalmiddlefrontal', 'lh-cuneus', 'lh-entorhinal', 'lh-fusiform', ...
    'lh-inferiorparietal', 'lh-inferiortemporal', 'lh-isthmuscingulate', ...
    'lh-lateraloccipital', 'lh-lateralorbitofrontal', 'lh-lingual', ...
    'lh-medialorbitofrontal', 'lh-middletemporal', 'lh-parahippocampal', ...
    'lh-paracentral', 'lh-parsopercularis', 'lh-parsorbitalis', ...
    'lh-parstriangularis', 'lh-pericalcarine', 'lh-postcentral', ...
    'lh-posteriorcingulate', 'lh-precentral', 'lh-precuneus', ...
    'lh-rostralanteriorcingulate', 'lh-rostralmiddlefrontal', ...
    'lh-superiorfrontal', 'lh-superiorparietal', 'lh-superiortemporal', ...
    'lh-supramarginal', 'lh-frontalpole', 'lh-temporalpole', ...
    'lh-transversetemporal', 'lh-insula', ...
    'rh-banks-of-superior-temporal-sulcus', 'rh-caudalanteriorcingulate', ...
    'rh-caudalmiddlefrontal', 'rh-cuneus', 'rh-entorhinal', 'rh-fusiform', ...
    'rh-inferiorparietal', 'rh-inferiortemporal', 'rh-isthmuscingulate', ...
    'rh-lateraloccipital', 'rh-lateralorbitofrontal', 'rh-lingual', ...
    'rh-medialorbitofrontal', 'rh-middletemporal', 'rh-parahippocampal', ...
    'rh-paracentral', 'rh-parsopercularis', 'rh-parsorbitalis', ...
    'rh-parstriangularis', 'rh-pericalcarine', 'rh-postcentral', ...
    'rh-posteriorcingulate', 'rh-precentral', 'rh-precuneus', ...
    'rh-rostralanteriorcingulate', 'rh-rostralmiddlefrontal', ...
    'rh-superiorfrontal', 'rh-superiorparietal', 'rh-superiortemporal', ...
    'rh-supramarginal', 'rh-frontalpole', 'rh-temporalpole', ...
    'rh-transversetemporal', 'rh-insula', 'left-cerebral-white-matter', ...
    'left-lateral-ventricle', 'left-inf-lat-vent', ...
    'left-cerebellum-white-matter', 'left-cerebellum-cortex', ...
    'left-thalamus-proper', 'left-caudate', 'left-putamen', ...
    'left-pallidum', 'third-ventricle', 'fourth-ventricle', 'brain-stem', ...
    'left-hippocampus', 'left-amygdala', 'csf', 'left-accumbens-area', ...
    'left-ventraldc', 'right-cerebral-white-matter', ...
    'right-lateral-ventricle', 'right-inf-lat-vent', ...
    'right-cerebellum-white-matter', 'right-cerebellum-cortex', ...
    'right-thalamus-proper', 'right-caudate', 'right-putamen', ...
    'right-pallidum', 'right-hippocampus', 'right-amygdala', ...
    'right-accumbens area', 'right-ventraldc'
};

region_map = containers.Map(all_regions, decoded_regions);

significant_result_count = 0;

% Initialize table to store significant results
struc_results = table('Size', [0 numel(varTypes)], 'VariableTypes', varTypes, 'VariableNames', colNames);

% Iterate through each predictor of interest
for pred_idx = 1:length(pred_vars)
    pred_struct = prop_stats.(pred_vars{pred_idx});

    % Iterate through each region
    for reg_idx = 1:numel(all_regions)

        region_name = all_regions{reg_idx};
        mdl_table = pred_struct.(region_name);

        new_row = {region_map(region_name), pred_vars{pred_idx}};
    
        % Iterate through each predictor
        for prop_idx = 1:numel(prop_names)
            prop_name = prop_names{prop_idx};
            mdl_prop_idx = find(ismember(mdl_table.Properties.RowNames, prop_name));

            % Add NaN values for subcortical thickness and white matter intensity (which don't exist)
            if (strcmp(prop_name, 'thick') || strcmp(prop_name, 'wm')) && reg_idx > 68
                new_row = [new_row, repmat({NaN}, 1, (numel(colNames) - 2) / 3)];

            % Check if beta, model, and intercept p-values are significant (p <= 0.05)
            elseif mdl_table.model_pvalue(mdl_prop_idx) <= 0.05 && mdl_table.qvalue(mdl_prop_idx) <= 0.05 && mdl_table.Intercept_pvalue(mdl_prop_idx) <= 0.05
                significant_result_count = significant_result_count + 1;

                % Create string of confidence intervals
                std_beta_CI = ['[', num2str(round(mdl_table.std_beta_CI_low(mdl_prop_idx), 4)), ...
                               ', ', num2str(round(mdl_table.std_beta_CI_up(mdl_prop_idx), 4)), ']'];
                beta_CI = ['[', num2str(round(mdl_table.beta_CI_low(mdl_prop_idx), 4)), ', ', ...
                           num2str(round(mdl_table.beta_CI_up(mdl_prop_idx), 4)), ']'];

                if strcmp(std_beta_CI, '[NaN, NaN]')
                    new_row = [new_row {mdl_table.regression_coef(mdl_prop_idx), beta_CI, mdl_table.qvalue(mdl_prop_idx)}];
                else
                    new_row = [new_row {mdl_table.std_beta(mdl_prop_idx), std_beta_CI, mdl_table.qvalue(mdl_prop_idx)}];
                end

                adj_tbl = all_stats.(pred_vars{pred_idx}).(all_regions{reg_idx}).(prop_names{prop_idx});
                adj_stats = cell(1, 4);
                
                for adj_idx = 1:numel(adj_tbl.Properties.RowNames)
                    % Get adjustment beta if significant
                    if adj_tbl.qvalue(adj_idx) <= 0.05 && ~contains(adj_tbl.Properties.RowNames{adj_idx}, '_x_')

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

                % Add interaction stats if significant
                if do_interaction
                    if adj_tbl.qvalue(interaction_idx) <= 0.05
                        interaction_stats = cell(1, 3);
                        interaction_stats{1} = num2str(round(adj_tbl.std_beta(interaction_idx), 4));
                        interaction_stats{2} = ['[', num2str(round(adj_tbl.std_beta_CI_low(interaction_idx), 4)), ...
                            ', ', num2str(round(adj_tbl.std_beta_CI_up(interaction_idx), 4)), ']'];
                        interaction_stats{3} = num2str(round(adj_tbl.qvalue(interaction_idx), 4));
                        new_row = [new_row interaction_stats];
                    else
                        new_row = [new_row, repmat({NaN}, 1, 3)];
                    end
                end

            else
                new_row = [new_row, repmat({NaN}, 1, (numel(colNames) - 2) / 3)];
            end
        end

        % Add row to significant results table
        if sum(cellfun(@(x) isnumeric(x) && isnan(x), new_row)) < numel(colNames) - 2
            struc_results = [struc_results; new_row];
        end
    end
end

% Save Excel file
output_dir = fileparts(outputfile);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
writetable(struc_results, outputfile, 'Sheet', tabname);
disp([num2str(significant_result_count), ' significant results.']);
disp(['Saved to ' outputfile]);

end