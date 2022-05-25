# Functions for Single-Case Methodologies
Matlab, Python and R functions for single-case methodologies


## Motivation
Various methodologies for singe-case studies have been developed by Dr. John R. Crawford and his colleagues (Crawford & Garthwaite, 2002; Crawford, Garthwaite & Porter, 2010; Crawford & Howell, 1998). 
While they provided computer programs for implementing these methodologies, these programs are executable files (.exe) for Micorsoft Windows OS and the source codes are not available.
Researchers using other operating systems (e.g., MacOS or Linux) cannot run these programs and implement these methods

Further, the programs require users to manually input values to run a single test. 
In other words, users need to calculate the required values (e.g., mean and standard deviation of the controls' scores) in other software (e.g., Excel, R, Python) and then enter them into the programs.
This is an error-prone process (i.e., may incorrectly enter the values) and is not feasible when users need to run analyses on multiple tests/measures

Therefore, I wrote functions of the singe-case methodologies in Matlab, Python and R based on the formulae in the papers by Crawford and his colleagues.

Relative to the Windows exe programs provided by Dr. John R. Crawford, there are several key strenghs of the functions here:
* Researchers can implement the single-case methodologies on computers of all kinds of OS, including Microsoft Windows, MacOS, Linux, as long as they have access to Matlab, Python, or R.
* Researchers can calculate the required values and perform single-case analyses in a single script
* Researchers can integrate the single-case methodologies into their analytical pipelines and can easily apply them to multiple tests/measures
* It is less error-prone since researchers do not need to manually enter values into the programs to perform the analysis


## Overview of the Functions
