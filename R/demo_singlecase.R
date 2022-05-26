# Demonstrate How to Use Crawford's Single-Case Methodologies in R

#### Preparation ####
# Clear workspace
rm(list=ls(all=TRUE))

# Load Necessary Package
library(tidyverse) 
library(pander) # for pander function: prints an R object in Pandoc's markdown

# Load functions for single-case methods
cwd <- getwd()
source(paste(cwd,"CrawfordHowell.R",sep="/"))
source(paste(cwd,"CrawfordRSDT.R",sep="/"))


#### Load Sample Data ####
# Load sample data of a patient and controls
filePath <- paste(dirname(cwd), "data_singlecase.csv", sep="/")
data_all <- read_csv(filePath)

# Separate patient's data from controls'
data_p <- data_all %>% filter(Group == "Patient") # the patient's data
data_c <- data_all %>% filter(Group == "Control") # controls' data


#### Single-Case T-Test ####
# Testing whether or not the patient's accuracy is significantly different from 
# the control group
Conds <- distinct(data_all, Condition)

# Conduct the single-case t-test for each condition separately
res_singleT <- tibble(
  tVal = numeric(),
  df = numeric(),
  pVal = numeric(),
  Zcc = numeric(),
  perct_ctrlBelowCase = numeric()
) # empty tibble/dataframe storing the results
for (i in 1:dim(Conds)[1]) {
  caseData <- data_p %>% 
    filter(Condition == pull(Conds[i,])) %>% 
    pull(Accuracy) # accuracy of the patient
  controlData <- data_c %>% 
    filter(Condition == pull(Conds[i,])) %>% 
    pull(Accuracy) # accuracies of controls
  res_singleT[i,] <- CrawfordHowell(caseData, controlData)
}
res_singleT <- res_singleT %>% 
  add_column(Condition = pull(Conds), .before = "tVal")

# Display the results
# cat("Crawford's T-test:")
res_singleT %>%
  mutate(across(where(is.numeric), round, digits=3)) %>%
  pander()


#### Revised standardized dissociation test (RSDT) ####
# Compare the difference between the patient's performance on two conditions of 
# the same task with the distribution of differences in controls.
# Different-view vs. Same-view conditions

# Patient's accuracies of the 2 conditions 
acc_p <- data_p %>% 
  pivot_wider(id_cols = c(Subject, Group), names_from = Condition, values_from = Accuracy) %>% 
  select(Diff, Same)
# Controls' accuracyies of the 2 conditions
acc_c <- data_c %>%
  pivot_wider(id_cols = c(Subject, Group), names_from = Condition, values_from = Accuracy) %>% 
  select(Diff, Same)

# RSDT
# cat("Crawford's Revised Standardized Dissociation Test (Different-view vs. Same-view):")
CrawfordRSDT(acc_p, acc_c) %>% 
  pander()

