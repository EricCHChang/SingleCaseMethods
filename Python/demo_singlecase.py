# %% [markdown]
# # Demonstrate How to Use Crawford's Single-Case Methodologies in Python

# %%
# Import Python packages
import pandas as pd
from pathlib import Path

# Import all functions for single-case methods
from crawfordtests.singlecase_t import crawford_t
from crawfordtests.singlecase_dissocs import crawford_udt, crawford_rsdt

# %% [markdown]
### Load Sample Data
#Load sample data of a patient and controls
# 
#The csv file contains data of a face oddity task from a single patient
#and his demographically-matched controls.
#
#There are two conditions of the face oddity task: different- and same-view conditions.
#
#The dependent variable/measure is accuracy (proportion of correct trials).
# 
#For details of the face oddity task, see Lee, A. C. H., Buckley, M. J., Pegman, S. J., Spiers, H., Scahill, V. L., Gaffan, D., ... & Graham, K. S. (2005). Specialization in the medial temporal lobe for processing of objects and scenes. *Hippocampus, 15*(6), 782-797.
#
#(note: the sample data here does not come from Lee et al., 2005 but from a different study.) 


# %%
# Load sample data of a patient and controls
cwd = Path.cwd()
filePath = cwd.parent / 'data_singlecase.csv'
data_all = pd.read_csv(filePath)
print(data_all.info())
data_all.head()

# %%
# Separate patient's data from controls'
data_p = data_all.loc[data_all['Group']=='Patient'] # the patient's data
data_c = data_all.loc[data_all['Group']=='Control'] # all controls' data

# %% [markdown]
# ## Single-Case T-Test
# Testing whether or not the patient's accuracy is significantly different from the control group

# %%
_res_df = pd.DataFrame()
# Conduct the single-case t-test for each condition separately
for cond in data_all['Condition'].unique().tolist():
    _tmp = crawford_t(
        case=data_p[data_p['Condition']==cond],
        control=data_c[data_c['Condition']==cond],
        dv='Accuracy'
    )
    _tmp.index = [cond]
    _res_df = pd.concat([_res_df, _tmp], axis=0)
    # print('Crawford T-test of condition', cond+':')
    # print(_tmp, '\n')

print("Crawford's T-test:")
_res_df.round(4)

# %% [markdown]
#**Summary**
#* The patient's accuracy is significantly lower than the control group in the different-view condition
#* In contrast, there is no difference in accuracy of the same-view condition between the patient and controls

# %% [markdown]
### Dissociation Test
#Compare the difference between the patient's performance on two conditions of the same task
#with the distribution of differences in controls.
# 
#* Unstandardized dissociation test (UDT)
#* Revised standardized dissociation test (RSDT)

# %%
# UDT
print("Crawford's Unstandardized Dissociation Test (Different-view vs. Same-view):")
crawford_udt(
    case=data_p, control=data_c, 
    dv='Accuracy', subject='Subject',
    cond='Condition', orderCond=['Diff','Same']
).round(4)

# %%
# RSDT
print("Crawford's Revised Standardized Dissociation Test (Different-view vs. Same-view):")
crawford_rsdt(
    case=data_p, control=data_c, 
    dv='Accuracy', subject='Subject',
    task='Condition', orderTask=['Diff','Same']
).round(4)

# %% [markdown]
#**Summary**
#* The patient's difference in accuracy between the two conditions is different from that of controls, but it is only significant when the scores are not standardized.

# %%
