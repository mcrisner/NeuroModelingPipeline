function compare_results(best_file, second_file, output_file, level_of_analysis)

% Description:
%   This function compares two Excel files containing model results by identifying rows 
%   from the first file ('best_file') that are also present in the second file ('second_file') 
%   across all shared sheets. The intersected rows are saved to a new Excel file.
%
% Input Arguments:
%   best_file         - String specifying the path to the first Excel file.
%   second_file       - String specifying the path to the second Excel file.
%   output_file       - String specifying the name of the output Excel file to write results to.
%   level_of_analysis - String indicating the level of analysis ('node', 'network', 'connectome') 
%                       to determine relevant columns for comparison.
%
% Process:
%   1. Retrieves all sheet names from both Excel files.
%   2. Identifies sheets common to both files.
%   3. Iterates over each shared sheet:
%      - Reads corresponding tables from both files.
%      - Subsets relevant columns based on the level of analysis.
%      - Identifies rows in the best_file that also appear in the second_file.
%      - If no intersection is found or the tables are empty, creates a blank sheet.
%   4. Writes the intersected (or blank) table to a new Excel file, preserving sheet names.
%
% Notes:
%   - The column sets compared depend on the specified level of analysis:
%       'node'      → {'Outcome', 'Predictor', 'Node'}
%       'network'   → {'Outcome', 'Predictor', 'Network'}
%       'connectome'→ {'Outcome', 'Predictor'}
%   - Assumes variable names in the input tables are preserved as-is.
%   - Outputs will only include rows from 'best_file' that are exactly matched in 'second_file'.

try
    % Get sheet names from each file
    sheets1 = sheetnames(best_file);
    sheets2 = sheetnames(second_file);
    
    % Ensure both files have the same sheets
    common_sheets = intersect(sheets1, sheets2);
    
    % Loop over each common sheet and perform the intersection
    for i = 1:length(common_sheets)
        sheet_name = common_sheets{i};
        
        % Read tables from each sheet
        best_results = readtable(best_file, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');
        second_results = readtable(second_file, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');
        
        if ~isempty(best_results) && ~isempty(second_results)
    
            % Select relevant columns for comparison
            if strcmp(level_of_analysis, 'node')
                best_subset = best_results(:, {'Outcome', 'Predictor', 'Node'});
                second_subset = second_results(:, {'Outcome', 'Predictor', 'Node'});
            elseif strcmp(level_of_analysis, 'network')
                best_subset = best_results(:, {'Outcome', 'Predictor', 'Network'});
                second_subset = second_results(:, {'Outcome', 'Predictor', 'Network'});
            elseif strcmp(level_of_analysis, 'connectome')
                best_subset = best_results(:, {'Outcome', 'Predictor'});
                second_subset = second_results(:, {'Outcome', 'Predictor'});
            end
            
            % Find the intersection
            [~, idx] = ismember(best_subset, second_subset, 'rows');
            intersection_rows = idx > 0;
            filtered_best_results = best_results(intersection_rows, :);
    
        else
    
            % Create a blank tab except for the column names
            filtered_best_results = array2table(NaN(0, width(best_results)), 'VariableNames', best_results.Properties.VariableNames);
        
        end
    
        % Write the results to a new file in the corresponding sheet
        writetable(filtered_best_results, output_file, 'Sheet', sheet_name);

    end

catch ME
    % Display error information
    fprintf('Insufficient results:\n');
    fprintf('%s\n', ME.message);
end

end
