import numpy as np
import pandas as pd
from math import sqrt
from scipy import stats
import scipy.special as special
# import warnings

def crawford_t(case, control, dv=None, alphaLevel=0.05, tail='two'):
    """Crawford-Howell T-test For Single Case Analysis 
    Compare the score of a single case to the matched control group
    
    See Crawford & Howell (1998) for the formula of the modified T-test 
    
    See Crawford & Garthwaite (2002) for the confidence limit of the percentage estimate 
    
    See Crawford, Garthwaite & Porter (2010) for the formula of the effect size 

    Args:
    ----
        case (pd.DataFrame or pd.Series): Data of the single case 
                                          (must be in wide format if it's a dataframe)
        
        control (pd.DataFrame or pd.Series): Data of the controls that are matched to the single case 
                                             (must be in wide format if it's a dataframe)
        
        dv (str): The name of the column containing the dependent variable being compared
                  (None if the case and control are pd.Series)
        
        alphaLevel (float, optional): Significance level (type 1 error). Default is 0.05.
        
        tail (str, optional): One-sided (`'one'`) or Two-sided (`'two'`) t-test. Defaults is `'two'`.

    Raises:
    ------
        ValueError: The value of ``dv`` is missing in the single case data.

    Returns:
    -------
        results_df (pd.Dataframe): Statistics of the Crawford-Howell T-test
            * ``'tVal'``: t-value
            * ``'df'``: degrees of freedom (number of controls - 1)
            * ``'pVal'``: p-value
            * ``'Zcc'``: effect size (the magnitude of the difference between the single case and the controls)
            * ``'Zcc_95CI'``: confidence interval (default 95%) of the effect size [lower bound, upper bound]
            * ``'normBelowCase%'``: percentage of normal/control population haveing a score below the case score
            * ``'normBelowCase%_95CI'``: confidence limit (default 95%) of the percentage estimate [lower bound, upper bound]
        
    """
        
    if dv is not None:
        # Check that case is a dataframe 
        assert isinstance(case, pd.DataFrame), 'Data of the single case must be a pandas dataframe when dv is provided.'
        # Check that control is a dataframe 
        assert isinstance(control, pd.DataFrame), 'Data of the controls must be a pandas dataframe when dv is provided.'
        # The dv must be a string
        assert isinstance(dv, str), 'dv must be a string'
    else:
        # Check that case is a pandas Series
        assert isinstance(case, pd.Series), 'Data of the single case must be a pandas series when dv is not provided.'
        # Check that case is a pandas Series
        assert isinstance(control, pd.Series), 'Data of the controls must be a pandas series when dv is not provided.'

        # dv is not provided and both case and control data are pd.Series
        # convert the case data to a dataframe
        case  = pd.melt(pd.DataFrame(case), value_name='dv')
        # convert the control data to a dataframe
        control  = pd.melt(pd.DataFrame(control), value_name='dv')

        dv = 'dv'
        
    # Check that case is has only 1 single subject 
    assert len(case)==1, 'There must be only 1 single case subject in the case dataframe.'
    
    # Check that the value of the dv exists in the case dataframe
    if case[dv].isnull().any():
        raise ValueError('The value of the variable ' + dv + ' is missing in the single case data.')
    
    # Remove missing values (NaN) in the controls data
    if control[dv].isnull().any():
        control = control.dropna(subset=[dv])
        # idx = control[dv].isnull()
        # control = control.loc[~idx]
        print('Missing values in the control. They have been removed.')
    
    # Calculate values needed for the test    
    nCtrl = len(control) # number of controls
    idx_case = case[dv].index[0] # in case the index is not zero
    case_val = case.loc[idx_case, dv] # value of the DV of the single case
    ctrl_avg = control[dv].mean() # mean of the DV in controls
    ctrl_sd = control[dv].std() # standard deviation of the DV in controls
    
    if ctrl_sd==0:
        print("Standard deviation of controls is zero. Return NaN")
        tVal = np.nan
        dof = nCtrl - 1
        pVal = np.nan
        effSize = np.nan
        effSize_CIlower = np.nan
        effSize_CIupper = np.nan
        percent_ctrlBelowCase = np.nan
        percent_ctrlBelowCase_CIlower = np.nan
        percent_ctrlBelowCase_CIupper = np.nan
    
    else:    
            
        # Crawford-Howell single case t-test
        tVal = (case_val - ctrl_avg) / ( ctrl_sd * sqrt((nCtrl+1)/nCtrl) ) # t-value
        dof = nCtrl - 1 # degrees of freedom
        pVal = stats.t.sf(abs(tVal), df=dof) * 2 # two-tailed
        # pVal = (1-stats.t.cdf(abs(tVal), df=dof)) * 2 # two-tailed
        if tail=='one':
            # report 1-tailed p-value instead
            pVal = pVal/2 # one-tailed
            
        # Effect size (as a z-score)
        effSize = (case_val - ctrl_avg) / ctrl_sd
        
        # Calculate the confidence interval of the effect size (as z-score)
        # see Appendix A in Crawford, Garthwaite & Porter (2010)
        nctVal = effSize * sqrt(nCtrl)
        nc_lower = special.nctdtrinc(dof, 1-(alphaLevel/2), nctVal)
        nc_upper = special.nctdtrinc(dof, alphaLevel/2, nctVal)
        # lower bound: 0.975th percentile non-centrality parameter divided by sqrt(nCtrl)
        # upper bound: 0.025th percentile non-centrality parameter divided by sqrt(nCtrl)
        effSize_CIlower = nc_lower / sqrt(nCtrl)
        effSize_CIupper = nc_upper / sqrt(nCtrl)

        # Calculate the estimated percentage of normal population falling below case's score
        percent_ctrlBelowCase = stats.norm.cdf(effSize) * 100 # multiply by 100 so that it's percentage

        # Calculate confidence limit of the percentage estimate 
        # (percentage of the control population exhibiting a lower score than the single case)
        # see section 2 and Fig.1 in Crawford & Garthwaite (2002)
        # *************** Values needed in this chunk have been calculated above ***************************
        # c1 = (case_val - ctrl_avg) / ctrl_sd # this is equal to effect size (Z-cc)
        # nctVal = c1 * sqrt(nCtrl) # t-value of a non-central t-distribution (Z-cc * sqrt(nCtrl))
        # # calculate the non-centrality parameter for non-central t distribution at given percentile
        # # scipy.special.nctdtrinc(df, p, t) # p: percentile; t: t-value
        # nc_lower = special.nctdtrinc(dof, 1-(alphaLevel/2), nctVal) # Crawford used delta in his paper
        # nc_upper = special.nctdtrinc(dof, alphaLevel/2, nctVal) 
        # ***************************************************************************************************
        # lower limit: probability of the point estimate of the percentage lower than (nc_lower/sqrt(nCtrl)) in a normal distribution
        # upper limit: probability of the point estimate of the percentage lower than (nc_upper/sqrt(nCtrl)) in a normal distribution
        percent_ctrlBelowCase_CIlower = stats.norm.cdf(nc_lower / sqrt(15)) * 100 # multiply by 100 so that it's percentage
        percent_ctrlBelowCase_CIupper = stats.norm.cdf(nc_upper / sqrt(15)) * 100
            
    # Return
    results_df = pd.DataFrame(
        columns=['tVal','df','pVal','Zcc','Zcc_95CI','normBelowCase%','normBelowCase%_95CI']
    )
    results_df.loc[0,:] = np.array(
        [tVal, dof, pVal,
        effSize,
        np.around([effSize_CIlower, effSize_CIupper], 4),
        percent_ctrlBelowCase,
        np.around([percent_ctrlBelowCase_CIlower, percent_ctrlBelowCase_CIupper], 4)],
        dtype='object'
    )
    results_df = results_df.astype(
        {'tVal':'float64', 'df':'int64', 'pVal':'float64',
         'Zcc':'float64', 'normBelowCase%':'float64'}
    )
    results_df.index = ['Ttest_cc']
    
    # results = {'tVal': tVal, 'df': dof, 'pVal': pVal, 'Zcc': effSize}
    # results_df = pd.DataFrame(results, index=['Ttest_cc'])

    return results_df
    
    
    # Other functions related to non-central t-distribution
    # Given a quntile and a non-centrality parameter, calculate the t-value for non-central t distribution
    # scipy.stats.nct.ppf(q, df, nc) # q: quantile; nc: non-centrality parameter

    # Given a non-centrality parameter, calculate the percentile for for non-central t distribution
    # scipy.special.nctdtr(df, nc, t) # nc: non-centrality parameter; t: t-value
        

    
    
    