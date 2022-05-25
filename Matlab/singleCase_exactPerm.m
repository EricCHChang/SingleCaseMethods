function [meanDiff_act,meanDiff_perm,rnk,p_val_Perm,CI_perm] = singleCase_exactPerm(singleCase,controlGroup,alphaLevel,tail)
% Exact permutation test for single-case studies
% 1. Compute the mean difference in dependent varialbe between the single 
% case (e.g., a patient) and the control group.
% 2. Out of all subjects, including the patient and controls, select one 
% control from the control group and treat it as the "patient".
% The real patient now is treated as control. Compute the mean difference.
% 3. Repeat step 2 for every control (i.e., draw without replacement).
% Eventually the number of mean differences will be equal to the number of 
% controls
% 4. See where the actual difference between the patient and controls fall
% in the distribution of all permuted mean differences. 
%
% Inputs:
%   singleCase - data from the single case subject (e.g., a patient)
%   controlGroup - data from the control group (a column vector)
%   alpha level - significance level (alpha)
%   tail - one-tailed (1 or -1) or two-tailed p-value (2)
%
% Outputs:
%   meanDiff_act - actual mean difference between the single case and the
%                  control group
%   meanDiff_perm - exact permuted mean difference. A 1 x nControls vector
%   rnk - rank of the actual mean difference among the permuted mean
%         difference
%   p_val - p value of the exact permutation test
%   CI_perm - confidence interval of the exact permutation test
%
% Written by C.-H. Eric Chang, Aug, 2019
%
%%

% Check missing data in the control group and if so, remove it
ind_missData = find(isnan(controlGroup));
nMissData = sum(isnan(controlGroup));
if nMissData==0
    nControls = size(controlGroup,1);
else
    nControls = size(controlGroup,1) - nMissData;
    controlGroup(ind_missData) = [];
    disp([num2str(nMissData) ' missing data in the control group'])
end

% Actual difference between the single case and the control group
meanDiff_act = singleCase - mean(controlGroup);
% meanDiff_act = mean(controlGroup) - singleCase;

% Permute all possible permutations
for i = 1:nControls
    singleCase_perm = controlGroup(i); % a permuted single case
    ind_controlRemained = setdiff(1:nControls, i);
    controlGroup_perm = [singleCase; controlGroup(ind_controlRemained)];
    meanDiff_perm(i) = singleCase_perm - mean(controlGroup_perm);  
%     meanDiff_perm(i) = mean(controlGroup_perm) - singleCase_perm; 
end

% Rank the actual mean difference among the permuted mean difference
if tail==2
    tmp_sort = sort(abs(meanDiff_perm), 'descend');
    tmp_act = abs(meanDiff_act);
    % how many permuted mean difference is greater or equal to the actual mean difference
    rnk = find(tmp_sort>=tmp_act, 1, 'last'); 
elseif tail==1
    tmp_sort = sort(meanDiff_perm, 'descend');
    tmp_act = meanDiff_act;
    % how many permuted mean difference is greater or equal to the actual mean difference
    rnk = find(tmp_sort>=tmp_act, 1, 'last');
elseif tail==-1
    tmp_sort = sort(meanDiff_perm, 'ascend');
    tmp_act = meanDiff_act;
    % how many permuted mean difference is greater or equal to the actual mean difference
    rnk = find(tmp_sort<=tmp_act, 1, 'last');
end
% the rank of the actual mean difference 
if isempty(rnk)
    rnk = 1;
else
    rnk = rnk + 1;
end

% Compute the p value
if tail==2 % default two-tailed test
    p_val_Perm = sum(abs(meanDiff_perm) >= abs(meanDiff_act)) / nControls;
elseif tail==1 % upper one-tailed test
    p_val_Perm = sum(meanDiff_perm >= meanDiff_act) / nControls;
elseif tail==-1 % lower one-tailed test
    p_val_Perm = sum(meanDiff_perm <= meanDiff_act) / nControls;
end

% Compute confidence interval of the permuted values
tmp_sort = sort(meanDiff_perm, 'descend');
CI_perm(1) = prctile(tmp_sort,100*alphaLevel/2); %CI lower bound
CI_perm(2) = prctile(tmp_sort,100-100*alphaLevel/2); % CI upper bound

end