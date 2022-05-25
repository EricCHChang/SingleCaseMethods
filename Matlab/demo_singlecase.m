% Demonstrate how to use Crawford's Single-Case Methodologies in Matlab

%% Load Sample Data
cwd = cd; % current directory
[parentPath, currentName, ext] = fileparts(cwd);
filePath = fullfile(parentPath, 'data_singlecase.csv');
data_all = readtable(filePath);

% Separate the patient's data from contros'
data_p = data_all(strcmp('Patient', data_all.Group), :); % patient's data
data_c = data_all(strcmp('Control', data_all.Group), :); % controls' data

% Names of the conditions
Conds = unique(data_all.Condition);

%% Single-Case T-Test
res_singleT = nan(length(Conds), 5); % empty matrix storing the results

% Conduct the single-case t-test for each condition separately
for i = 1:length(Conds)
    cond = Conds{i}; % name of the condition
    singleCase = data_p(strcmp(data_p.Condition, cond), :).Accuracy; % patient's accuracy
    controlGroup = data_c(strcmp(data_c.Condition, cond), :).Accuracy; % contros' accuracies
    % Crawford's t-test
    [tVal, df, pVal, effSize, perct_ctrlBelowCase] = CrawfordHowell(singleCase, controlGroup);
    % Store the results
    res_singleT(i, :) = [tVal, df, pVal, effSize, perct_ctrlBelowCase];
    
    clear tVal df pVal effSize perct_ctrlBelowCase
end

% Convert the result matrix to a table for readability
tbl_singleT = array2table(res_singleT, 'VariableNames',{'tVal','df','pVal','Zcc','perct_normBelowCase'});
tbl_singleT = addvars(tbl_singleT, Conds, 'Before','tVal');
tbl_singleT

%% Single-Case Exact Permutation Test
res_exactPerm = nan(length(Conds), 4); % empty matrix storing the results

% Conduct the single-case exact permutation test for each condition separately
for i = 1:length(Conds)
    cond = Conds{i}; % name of the condition
    singleCase = data_p(strcmp(data_p.Condition, cond), :).Accuracy; % patient's accuracy
    controlGroup = data_c(strcmp(data_c.Condition, cond), :).Accuracy; % contros' accuracies
    % Single-case exact permutation test
    [meanDiff_act,meanDiff_perm,rnk,p_val_Perm,CI_perm] = singleCase_exactPerm(singleCase, controlGroup);
    % Store the results
    res_exactPerm(i, :) = [meanDiff_act, p_val_Perm, CI_perm(1), CI_perm(2)];
    
    clear meanDiff_act meanDiff_perm rnk p_val_Perm CI_perm
end

% Convert the result matrix to a table for readability
tbl_exactPerm = array2table(res_exactPerm, 'VariableNames',{'meanDiff_act','pVal','CIofDiff_lower','CIofDiff_upper'});
tbl_exactPerm = addvars(tbl_exactPerm, Conds, 'Before','meanDiff_act');
tbl_exactPerm



