# Functions for Single-Case Methodologies
Python, R and Matlab functions for single-case methodologies



## Motivation
Various methodologies for singe-case studies have been developed by Dr. John R. Crawford and his colleagues (Crawford & Garthwaite, 2002, 2005; Crawford, Garthwaite & Porter, 2010; Crawford & Howell, 1998). 
While they provided computer programs for implementing these methodologies (see https://homepages.abdn.ac.uk/j.crawford/pages/dept/SingleCaseMethodsComputerPrograms.HTM), these programs are executable files (.exe) for Micorsoft Windows operating system (OS) and the source codes are not available.
Researchers using other operating systems (e.g., MacOS or Linux) cannot run these programs and implement these methods.

Further, the programs require users to manually input values to run a single test. 
In other words, users need to calculate the mean and standard deviation of controls' scores (and/or other statistics) in other software (e.g., Excel, R, Python) and then enter them into the programs.
This is an error-prone process (i.e., may incorrectly enter the values) and users will need to repeat the same steps multiple times for different tests/measures.

Therefore, I wrote functions of the singe-case methodologies in Matlab, Python and R based on the formulae in the papers by Crawford and his colleagues.

Relative to the Windows exe programs provided by Dr. John R. Crawford, key strength of the functions here include:
* Users can implement the single-case methodologies on computers of all kinds of OS, including Microsoft Windows, MacOS, Linux, as long as they have access to Matlab, Python, or R.
* Users can calculate the required values/statistics and perform single-case analyses in a single script.
* Users can integrate the single-case methodologies into their analytical pipelines and can easily apply them to multiple tests/measures.
* It is less error-prone since users do not need to manually enter values into the programs to perform analyses.



## Overview of the Functions

### Python

* `singlecase_t.py`: Crawford's t-test. This script corresponds to Crawford's `Singlims_ES.exe` program.

* `singlecase_dissocs.py`: Crawford's dissociation test, including the unstandardized difference test (UDT) and the revised standardized difference test (RSDT). This scirpt corresponds to Crawford's `Dissocs_ES.exe` and `RSDT_ES.exe` programs. Note that the python functions currently does not support calculating the interval estimate for the effect size and the a patient's abnormality of differences between two tests. 

### R 

* `CrawfordHowell.R`: Crawford's t-test. This script corresponds to Crawford's `Singlims_ES.exe` program. Note that this R function currently does not support calculating the interval estimate for the effect size and the abnormality of a single-case patient's score. 

* `CrawfordRSDT.R`: Crawford's revised standardized difference test (RSDT). This scirpt corresponds to Crawford's `RSDT_ES.exe` programs. Note that the functions currently does not support computing the interval estimate for the effect size and the a patient's abnormality of differences between two tests. 


### Matlab

* `CrawfordHowell.m`: Crawford's t-test. This script corresponds to Crawford's `Singlims_ES.exe` program. Note that this Matlab function currently does not support calculating the interval estimate for the effect size and the abnormality of a single-case patient's score.

* `singleCase_exactPerm.m`: A non-parametrics statistical test computing the significance of difference in a test (e.g., a neuropsychological test) between a single patient and a control sample via a exact permuation test. That is, it examines whether or not a patient exhibits a statistically significant deficit on a given test. Note that this is not part of Crawford's single-case methodologies. 



## Examples & Demonstration

### Python

* See the Jupyter Notebook [demo_singlecase.ipynb](Python/demo_singlecase.ipynb)



## Important Note

* The current functions for Crawford's single-case methodologies do not support the Bayesian approach (Crawford & Garthwaite, 2007).

* The results of the interval estimate computed by the Python functions are slightly different from those computed by Crawford's programs. It is likely due to differences in how the noncentrality parameter of a noncentral t-distribution is computed between different programming languages.



## References

* Crawford, J. R., & Garthwaite, P. H. (2002). Investigation of the single case in neuropsychology: Confidence limits on the abnormality of test scores and test score differences. *Neuropsychologia, 40*(8), 1196-1208.

* Crawford, J. R., & Garthwaite, P. H. (2005). Testing for suspected impairments and dissociations in single-case studies in neuropsychology: evaluation of alternatives using monte carlo simulations and revised tests for dissociations. *Neuropsychology, 19*(3), 318.

* Crawford, J. R., & Garthwaite, P. H. (2007). Comparison of a single case to a control or normative sample in neuropsychology: Development of a Bayesian approach. *Cognitive Neuropsychology, 24*(4), 343-372.

* Crawford, J. R., Garthwaite, P. H., & Porter, S. (2010). Point and interval estimates of effect sizes for the case-controls design in neuropsychology: rationale, methods, implementations, and proposed reporting standards. *Cognitive neuropsychology, 27*(3), 245-260.

* Crawford, J. R., & Howell, D. C. (1998). Comparing an individual's test score against norms derived from small samples. *The Clinical Neuropsychologist, 12*(4), 482-486.



