function mdl_stats = get_mdl_stats(mdl, ind, stat_lst, name_dim)
%   PURPOSE: get a set of common statistics for a model/predictor of interest
%
%   INPUTS:
%   mdl: linear or logistic model object
%   ind: index of predictor of interest in coefficient table (+1 to account
%   for intercept row)
%   stat_lst: optional list of statistics you would like returned
%   (otherwise, will default to return all statistics)
%   name_dim: optional: name row based on response variable (2, default), or
%   predictor of interest (1)
%
%   OUTPUTS:
%   mdl_stats: single table row of (requested) statistics for this model
%
%   NOTES: currently, subsetting by stat_lst is done at the end for
%   simplicity, since the stats are inexpensive. If adding expensive stats,
%   or running for large number of models, it may be worth it to avoid
%   running non-requested stats


%% setup variables and storage

%all available statistics
stats_all = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
    'std_beta', 'std_beta_CI_low', 'std_beta_CI_up','OR','OR_CI_low','OR_CI_up',...
    'SE', 'Wald', 'pvalue','Intercept_pvalue','model_pvalue','R_squared'};

%stats we will be returning:
if nargin < 3 || isempty(stat_lst)
    %most common stats as default
    stat_lst= {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
        'std_beta', 'std_beta_CI_low', 'std_beta_CI_up','OR','OR_CI_low','OR_CI_up',...
        'pvalue','Intercept_pvalue','model_pvalue','R_squared'};
else %in case stats we do not have or are not available are requested:
    if ~isempty(setdiff(stat_lst, stats_all))
        warning('some stats requested but cannot be calculated')
    end
    stat_lst = intersect(stat_lst, stats_all, 'stable'); 
end

%setup index of predictor of interest in variable table (may differ from
%coefficient table)
var_ind = find(strcmp(mdl.VariableNames, mdl.CoefficientNames(ind)));

%setup table of stats to be returned - row name as requested
if nargin < 4 || isempty(name_dim)
    name_dim = 2;
end
   
if name_dim==1
    row_lbl = mdl.VariableNames{var_ind}; 
elseif name_dim==2
    row_lbl = mdl.ResponseName;
end
mdl_stats= array2table(NaN(1, length(stat_lst)), ...
    'RowNames', {row_lbl}, 'VariableNames', stat_lst);

%% fill table
%get statistics from coefficient table
mdl_stats{1, 'regression_coef'} = mdl.Coefficients.Estimate(ind);
cis = mdl.coefCI;
mdl_stats{1, {'beta_CI_low', 'beta_CI_up'}} = cis(ind, :);
mdl_stats{1, 'SE'} = mdl.Coefficients.SE(ind);
tstat_val = mdl.Coefficients.tStat(ind);
mdl_stats{1, 'Wald'} = tstat_val^2;
mdl_stats{1, 'pvalue'} = mdl.Coefficients.pValue(ind);
mdl_stats{1, 'Intercept_pvalue'} = mdl.Coefficients.pValue(1);
mdl_stats{1,'model_pvalue'} = coefTest(mdl); 
mdl_stats{1,'R_squared'} = mdl.Rsquared.Adjusted;

%standardize beta value and confidence interval if both predictor of
%interest and response are not binary:
bin_resp = length(unique(rmmissing(mdl.Variables{:, mdl.ResponseName})))==2;
bin_pred = length(unique(rmmissing(mdl.Variables{:, var_ind})))==2;
if ~bin_resp && ~bin_pred
    st_dev_resp = std(mdl.Variables{:, mdl.ResponseName}, 'omitnan'); %standard deviation of response variable
    st_dev_pred = std(mdl.Variables{:, var_ind}, 'omitnan'); %standard deviation of predictor of interest

    st_beta = mdl.Coefficients.Estimate(ind)*(st_dev_pred/st_dev_resp);
    mdl_stats{1, 'std_beta'} = st_beta;

    std_cis= cis(ind,:) .* (st_dev_pred/st_dev_resp);
    mdl_stats{1,{'std_beta_CI_low', 'std_beta_CI_up'}} = std_cis; 
end

% give odds ratio values if binary response or predictor
if (bin_resp || bin_pred)
    mdl_stats{1, 'OR'} = exp(mdl_stats{1, 'regression_coef'});
    mdl_stats{1, 'OR_CI_low'} = exp(mdl_stats{1, 'beta_CI_low'});
    mdl_stats{1, 'OR_CI_up'} = exp(mdl_stats{1, 'beta_CI_up'});
end


%keep only what is requested
mdl_stats = mdl_stats(:, stat_lst);

end
