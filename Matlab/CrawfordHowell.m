function [tVal, df, pVal, effSize, perct_ctrlBelowCase] = CrawfordHowell(singleCase, controlGroup, alphaLevel, tail)
% Crawford-Howell t-test for single case studies
% 
% See Crawford & Howell (1998) for the formula of the modified T-test 
% See Crawford & Garthwaite (2002) for the confidence limit of the percentage estimate 
% See Crawford, Garthwaite & Porter (2010) for the formula of the effect size 
% 
% Inputs
%   singleCase - data from the single case subject (e.g., a single score of a patient)
%   controlGroup - data from the control group (a column vector)
% 
% Outputs
%   tVal - t-value of the 2-tailed test
%   df - degrees of freedom
%   pVal - p-value of the 2-tailed test
%
% Writtien by C.-H. Eric Chang, Oct 2017

%% Set default values
if ~exist('alphaLevel', 'var')
    alphaLevel = 0.95;
end

if ~exist('tail', 'var')
    tail = 'two';
end

%% Crawford-Howell single case t-test
% T-test
nContSubs = length(controlGroup);
tVal = (singleCase - mean(controlGroup)) ./ (std(controlGroup).*sqrt((nContSubs+1) ./ nContSubs));
df = nContSubs-1;
pVal = 2*(1-tcdf(abs(tVal),df));
if strcmp(tail, 'one')
    pVal = pVal/2;
end

%% Effect size
% Effect size (as a z-score)
effSize = (singleCase - mean(controlGroup)) / std(controlGroup);

%% Estimate percentage of the normal population falling below the single case's score
% Calculate the estimated percentage of normal population falling below case's score
mu = 0;
sigma = 1;
perct_ctrlBelowCase = cdf('Normal',effSize,mu,sigma) * 100; % multiply by 100 so that it's percentage

end