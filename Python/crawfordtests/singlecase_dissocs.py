# import numpy as np
import pandas as pd
import numpy as np
from math import sqrt
from scipy import stats

def crawford_udt(case, control, dv=None, cond=None, subject=None, 
                 orderCond=None, alphaLevel=0.05, tail='two'):
    """Crawford Dissociation Test For Single Case Analysis
    (a.k.a Unstandardized Difference Test, UDT)
    Compare the difference between a single case's performance on two conditions of the same task
    with the distribution of differences in controls.
    
    Note: the scores of the single case are not standardized (i.e., not converted to z-scores
    based on the control's performance) in this test, unlike the RSDT (see below).
    
    For example, when a single case (e.g., a patient) shows deficts in both condition A and B 
    (relative to controls), but the defict is larger for condition A than B, we may claim that 
    the single case has a strong dissociation between the two conditions. 
    
    See Crawford & Garthwaite (2005 for the formula of the UDT (equation 5)
    See Crawford, Garthwaite & Porter (2010) for the formula of the effect size of the difference

    Args:
        case (pd.DataFrame): Data of the single case (can be either long or wide format)
        
        control (pd.DataFrame): Data of the controls that are matched to the single case
                                (can be either long or wide format) 
        
        dv (str): The name of the column containing the dependent variable being compared. 
                  Default is None (only required when the data is in long format).
        
        cond (str, optional): The name of the column containing the conditions being compared. 
                              Default is None (only required when the data is in long format).
        
        subject (str, optional): The name of the column containing the subject IDs.
                                 Default is None (only required when the data is in long format).
        
        orderCond (list, optional): The list containing the order of the conditions being compared 
                                    (e.g., condition A - condition B or the opposite).
                                    Default is None (i.e., use the order of appearance of the conditoins 
                                    in the data).
        
        alphaLevel (float, optional): Significance level (type 1 error). Default is 0.05.
        
        tail (str, optional): One-sided or Two-sided t-test. Default is `'two'`.

    Raises:
        ValueError: If the data is long format, ``dv``, ``cond``, and ``subject`` arguments must be provided.
        ValueError: The value of ``dv`` is missing in the single case data.

    Returns:
        results_df (pd.Dataframe): Statistics of the Unstandardized Difference Test
            * ``'tVal'``: t-value
            * ``'df'``: degrees of freedom (number of controls - 1)
            * ``'pVal'``: p-value
            * ``'Zdcc'``: effect size (the magnitude of the difference between conditions in the single case
                          when it is compared to the difference in the controls)
    """

    
    # Check that case is a dataframe
    assert isinstance(case, pd.DataFrame), 'Data of the single case must be a pandas dataframe.'
    
    # Check that control is a dataframe
    assert isinstance(control, pd.DataFrame), 'Data of the controls must be a pandas dataframe.'
    
    # significance level must be between 0 and 1
    assert 0 < alphaLevel < 1, 'alphaLevel must be between 0 and 1.'
    
    # tail must be either one or two
    assert tail in ['two','one'], 'tail must be one or two.'
    
    # If the case and control are in wide format, convert them to long format
    if all([arg is None for arg in [dv, cond, subject]]):         
        # case
        # case_raw = case.copy()
        tmp_df = case.select_dtypes(include='number')
        assert tmp_df.shape[0]==1, 'There must be only 1 single case subject in the case dataframe.'
        assert tmp_df.shape[0]==2, 'There must be only 2 conditions being compared. Check the case data.'
        tmp_df['subject'] = np.arange(tmp_df.shape[0])
        case = pd.melt(tmp_df, id_vars='subject', var_name='cond', value_name='dv')
        del tmp_df

        # control
        # control_raw = control.copy()
        tmp_df = control.select_dtypes(include='number')
        assert tmp_df.shape[0] > 2, 'There must be at least 3 control subjects in the subject dataframe.'
        assert tmp_df.shape[0]==2, 'There must be only 2 conditions being compared. Check the control data.'
        tmp_df['subject'] = np.arange(tmp_df.shape[0])
        control = pd.melt(tmp_df, id_vars='subject', var_name='cond', value_name='dv')
        del tmp_df

        dv, cond, subject = 'dv', 'cond', 'subject'
    
    # One of the key arguments is missing when case and control are long format
    if any([arg is None for arg in [dv, cond, subject]]):
        raise ValueError('If the data is long format, dv, cond, and subject arguments must be provided.')
    
    # There must be only 2 different conditions to be compared
    assert len(case[cond].unique())==2, 'There must be only 2 conditions being compared. Check the case data.'
    assert len(control[cond].unique())==2, 'There must be only 2 conditions being compared. Check the control data.'
    
    # Check that case is has only 1 single subject (i.e., only 1 row for each condition)
    assert len(case)==2, 'There must be only 1 single case subject in the case dataframe.'
    
    # Check that the value of the dv exists in the case dataframe
    if case[dv].isnull().any():
        raise ValueError('There cannot be missing values in the variable ' + dv + ' in the single case data.')
    
    # Names of the 2 conditions 
    # (sorted in case the order of appearance is different in case and control)
    case_condNames = sorted(case[cond].unique())
    ctrl_condNames = sorted(control[cond].unique())
    assert case_condNames == ctrl_condNames, 'Names of conditions are not consistent between case and control.'
    condNames = sorted(case[cond].unique())    
    if orderCond: # if the order of the condition is specified
        assert isinstance(orderCond, list), 'orderCond must be a list'
        assert len(orderCond) == 2, 'There must be only 2 conditions in orderCond.'
        
        # check whether the condition names in the specified order are consistent with the data
        assert set(orderCond)==set(condNames), 'The condition names in orderCond are not consistent with those in the data.'
        # use the specified condition instead
        condNames = orderCond
    
    # Convert the control df to wide format
    control_wide = pd.pivot_table(control, values=dv, index=subject, columns=cond)
    
    # Remove missing values (NaN) in the control data 
    # pairwise removal: remove values in both conditions for a subjet even the value is only missing in 1 condition
    if control[dv].isnull().any():  
        # remvoe missing values (pairwise)
        control_wide = control_wide.dropna() #subset=[dv]
        # convert the control df back to long format
        control = pd.melt(control_wide, var_name=cond, value_name=dv, ignore_index=False).reset_index()
        print('Missing values in the control. They have been removed.')
    
    # Calculate values needed for the test 
    # (1) number of controls
    nCtrl = len(control[subject].unique())
    # (2) mean for each condition in the control
    ctrl_avg = control.groupby(by=cond)[dv].mean().reset_index()
    cond0_ctrl = ctrl_avg.loc[ctrl_avg[cond]==condNames[0], dv].values
    cond1_ctrl = ctrl_avg.loc[ctrl_avg[cond]==condNames[1], dv].values
    # (3) SD for each condition in the control
    ctrl_sd = control.groupby(by=cond)[dv].std().reset_index()
    cond0_sd_ctrl = ctrl_sd.loc[ctrl_sd[cond]==condNames[0], dv].values
    cond1_sd_ctrl = ctrl_sd.loc[ctrl_sd[cond]==condNames[1], dv].values
    # ctrl_var = control.groupby(by=cond)[dv].var().reset_index()
    # sum of variance of the 2 conditions in control
    sov = cond0_sd_ctrl**2 + cond1_sd_ctrl**2 #ctrl_var[dv].sum() 
    # (4) covariance of the 2 conditions in the control
    assert control_wide.cov().shape == (2,2), 'The covariance matrix should be 2-by-2.'
    ctrl_cov = control_wide.cov().iloc[0,1] 
    # (5) score of the single case in each condtion
    # the order of conditions must be identical to the control's conterpart
    cond0_case = case.loc[case[cond]==condNames[0], dv].values
    cond1_case = case.loc[case[cond]==condNames[1], dv].values    
    # (6) difference between the case and control, in each condition
    diff0 = cond0_case - cond0_ctrl    
    diff1 = cond1_case - cond1_ctrl
    
    # Crawford Unstandardized Difference Test     
    tVal = (diff0 - diff1) / sqrt( (sov-2*ctrl_cov) * ((nCtrl+1)/nCtrl) )
    dof = nCtrl - 1 # degrees of freedom
    pVal = stats.t.sf(abs(tVal), df=dof) * 2 # two-tailed
    # pVal = (1-stats.t.cdf(abs(tVal), df=dof)) * 2 # two-tailed
    if tail=='one':
        # report 1-tailed p-value instead
        pVal = pVal/2 # one-tailed
    
    # Effect size (as a z-score)
    if (cond0_sd_ctrl==0) | (cond1_sd_ctrl==0):
        print("Standard deviation of controls is zero in at least on task. Return effect size as NaN")
        effSize = np.nan
    else:
        corVal = control_wide.corr().iloc[0,1] # correlation between the 2 conditions in control
        effSize = ( (diff0/cond0_sd_ctrl) - (diff1/cond1_sd_ctrl) ) / sqrt(2-2*corVal)
    
    # Return
    results = {'tVal': tVal, 'df': dof, 'pVal': pVal, 'Zdcc': effSize}
    results_df = pd.DataFrame(results, index=['UDT'])

    return results_df
        
        
    
def crawford_rsdt(case, control, dv=None, task=None, subject=None,
                  orderTask=None, alphaLevel=0.05, tail='two'):
    """Crawford Revised Standardized Difference Test (RSDT) For Single Case Analysis
    Compare the difference between a single case's performance on two different tasks
    with the distribution of differences in controls.
    
    Note: the scores of the single case are standardized (i.e., converted to z-score based on
    the control's performance) in this test as the scale may be largely different in the two tasks.
    
    For example, when a single case (e.g., a patient) shows deficts in both task A and B 
    (relative to controls), but the defict is larger for task A than B, we may claim that 
    the single case has a strong dissociation between the two tasks. 
    
    See Crawford & Garthwaite (2005 for the formula of the RSDT (equations 6, 7, 8, 9)
    See Crawford, Garthwaite & Porter (2010) for the formula of the effect size of the difference

    Args:
        case (pd.DataFrame): Data of the single case (can be either long or wide format)
        
        control (pd.DataFrame): Data of the controls that are matched to the single case
                                (can be either long or wide format) 
        
        dv (str): The name of the column containing the dependent variable being compared. 
                  Default is None (only required when the data is in long format).
                  
        task (str, optional): The name of the column containing the tasks being compared. 
                              Default is None (only required when the data is in long format).
                              
        subject (str, optional): The name of the column containing the subject IDs.
                                 Default is None (only required when the data is in long format).
        
        orderTask (list, optional): The list containing the order of the tasks being compared 
                                    (e.g., task A - task B or the opposite).
                                    Default is None (i.e., use the order of appearance of the tasks
                                    in the data).
                                    
        alphaLevel (float, optional): Significance level (type 1 error). Default is 0.05.
        
        tail (str, optional): One-sided or Two-sided t-test. Default is `'two'`.

    Raises:
        ValueError: If the data is long format, ``dv``, ``task``, and ``subject`` arguments must be provided.
        ValueError: The value of ``dv`` is missing in the single case data.

    Returns:
        results_df (pd.Dataframe): Statistics of the Revised Standardized Difference Test
            * ``'tVal'``: t-value
            * ``'df'``: degrees of freedom (number of controls - 1)
            * ``'pVal'``: p-value
            * ``'Zdcc'``: effect size (the magnitude of the difference between tasks in the single case
                          when it is compared to the difference in the controls)
    """
    
    
    # Check that case is a dataframe
    assert isinstance(case, pd.DataFrame), 'Data of the single case must be a pandas dataframe.'
    
    # Check that control is a dataframe
    assert isinstance(control, pd.DataFrame), 'Data of the controls must be a pandas dataframe.'
    
    # significance level must be between 0 and 1
    assert 0 < alphaLevel < 1, 'alphaLevel must be between 0 and 1.'
    
    # tail must be either one or two
    assert tail in ['two','one'], 'tail must be one or two.'
    
    # If the case and control are in wide format, convert them to long format
    if all([arg is None for arg in [dv, task, subject]]):         
        # case
        # case_raw = case.copy()
        tmp_df = case.select_dtypes(include='number')
        assert tmp_df.shape[0]==1, 'There must be only 1 single case subject in the case dataframe.'
        assert tmp_df.shape[0]==2, 'There must be only 2 tasks being compared. Check the case data.'
        tmp_df['subject'] = np.arange(tmp_df.shape[0])
        case = pd.melt(tmp_df, id_vars='subject', var_name='task', value_name='dv')
        del tmp_df

        # control
        # control_raw = control.copy()
        tmp_df = control.select_dtypes(include='number')
        assert tmp_df.shape[0] > 2, 'There must be at least 3 control subjects in the subject dataframe.'
        assert tmp_df.shape[0]==2, 'There must be only 2 tasks being compared. Check the control data.'
        tmp_df['subject'] = np.arange(tmp_df.shape[0])
        control = pd.melt(tmp_df, id_vars='subject', var_name='task', value_name='dv')
        del tmp_df

        dv, task, subject = 'dv', 'task', 'subject'
   
    # One of the key arguments is missing when case and control are long format
    if any([arg is None for arg in [dv, task, subject]]):
        raise ValueError('If the data is long format, dv, task, and subject arguments must be provided.')
        
    # There must be only 2 different tasks to be compared
    assert len(case[task].unique())==2, 'There must be only 2 tasks being compared. Check the case data.'
    assert len(control[task].unique())==2, 'There must be only 2 tasks being compared. Check the control data.'

    # Check that case is has only 1 single subject (i.e., only 1 row for each condition)
    assert len(case)==2, 'There must be only 1 single case subject in the case dataframe.'

    # Check that the value of the dv exists in the case dataframe
    if case[dv].isnull().any():
        raise ValueError('There cannot be missing values in the variable ' + dv + ' in the single case data.')
    
    # Names of the 2 tasks
    # (sorted in case the order of appearance is different in case and control)
    case_taskNames = sorted(case[task].unique())
    ctrl_taskNames = sorted(control[task].unique())
    assert case_taskNames == ctrl_taskNames, 'Names of tasks are not consistent between case and control.'
    taskNames = sorted(case[task].unique())
    if orderTask: # if the order of the task is specified
        assert isinstance(orderTask, list), 'orderTask must be a list'
        assert len(orderTask) == 2, 'There must be only 2 tasks in orderTask.'
        
        # check whether the condition names in the specified order are consistent with the data
        assert set(orderTask)==set(taskNames), 'The task names in orderTask are not consistent with those in the data.'
        # use the specified condition instead
        taskNames = orderTask
    
    # Convert the control df to wide format
    control_wide = pd.pivot_table(control, values=dv, index=subject, columns=task)

    # Remove missing values (NaN) in the control data 
    # pairwise removal: remove values in both tasks for a subject even the value is only missing in 1 task
    if control[dv].isnull().any():  
        # remvoe missing values (pairwise)
        control_wide = control_wide.dropna() #subset=[dv] 
        # convert the control df back to long format
        control = pd.melt(control_wide, var_name=task, value_name=dv, ignore_index=False).reset_index()
        print('Missing values in the control. They have been removed.')
    
    # Calculate values needed for the test 
    # (1) number of controlspai
    nCtrl = len(control[subject].unique())
    # (2) mean for each condition in the control
    ctrl_avg = control.groupby(by=task)[dv].mean().reset_index()
    cond0_ctrl = ctrl_avg.loc[ctrl_avg[task]==taskNames[0], dv].values
    cond1_ctrl = ctrl_avg.loc[ctrl_avg[task]==taskNames[1], dv].values
    # (3) SD for each condition in the control
    ctrl_sd = control.groupby(by=task)[dv].std().reset_index()
    cond0_sd_ctrl = ctrl_sd.loc[ctrl_sd[task]==taskNames[0], dv].values
    cond1_sd_ctrl = ctrl_sd.loc[ctrl_sd[task]==taskNames[1], dv].values
    # (4) correlation between the 2 tasks in control
    corVal = control_wide.corr().iloc[0,1]
    # (5) score of the single case in each condtion
    # the order of tasks must be identical to the control's conterpart
    cond0_case = case.loc[case[task]==taskNames[0], dv].values
    cond1_case = case.loc[case[task]==taskNames[1], dv].values 
    # (6) difference between the case and control, in each task
    diff0 = cond0_case - cond0_ctrl
    diff1 = cond1_case - cond1_ctrl
    # (7) critical value (two-tailed) for t on n-1 degrees of freedom
    dof = nCtrl - 1 # degrees of freedom
    critVal = stats.t.isf(alphaLevel/2, df=dof)
    # (8) parameters needed to calculate the t-value of the RSDT
    if (cond0_sd_ctrl==0) | (cond1_sd_ctrl==0):
        print("Standard deviation of controls is zero in at least on task. Return NaN")
        tVal = np.nan
        dof = nCtrl - 1
        pVal = np.nan
        effSize = np.nan
        
    else:
        a = (1+corVal)*(1-corVal**2)
        b = (1-corVal)*( 4*(nCtrl-1)**2 + 4*(1+corVal)*(nCtrl-1) + (1+corVal)*(5+corVal) )
        c = -2 * ( (diff0/cond0_sd_ctrl) - (diff1/cond1_sd_ctrl) )**2 * ( (nCtrl*(nCtrl-1)**2)/(nCtrl+1) )
        
        # Crawford Revised Standardized Difference Test (RSDT)
        denominator = ((nCtrl+1)/nCtrl) * ( (2-2*corVal) + (2*(1-corVal**2))/(nCtrl-1) + ((5+critVal**2)*(1-corVal**2))/(2*(nCtrl-1)**2) +\
                                            (corVal*(1+critVal**2)*(1-corVal**2))/(2*(nCtrl-1)**2) )
        denominator = sqrt(denominator)
        psi = ( (diff0/cond0_sd_ctrl) - (diff1/cond1_sd_ctrl) ) / denominator
        
        tVal = sqrt( (-b + sqrt(b**2 - 4*a*c)) / (2*a) )
        tVal = tVal if psi>0 else -tVal # add sign (direction) to the t-value based on the psi value
        pVal = stats.t.sf(abs(tVal), df=19) * 2 # two-tailed
        # pVal = (1-stats.t.cdf(abs(tVal), df=dof)) * 2 # two-tailed
        if tail=='one':
            # report 1-tailed p-value instead
            pVal = pVal/2 # one-tailed
            
        # h = True if abs(psi)>critVal else False # whether significant or not

        # Effect size (as a z-score)
        effSize = ( (diff0/cond0_sd_ctrl) - (diff1/cond1_sd_ctrl) ) / sqrt(2-2*corVal)

    # Return
    results = {'tVal': tVal, 'df': dof, 'pVal': pVal, 'Zdcc': effSize}
    results_df = pd.DataFrame(results, index=['RSDT'])

    return results_df