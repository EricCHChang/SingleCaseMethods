#' Crawford's Revised Standardized Difference Test (RSDT) for testing differences 
#' among a patient's scores on two different tasks
#' The formula is based on Crawford, J. R., & Garthwaite, P. H. (2005). 
#' Testing for suspected impairments and dissociations in single-case studies in 
#' neuropsychology: evaluation of alternatives using monte carlo simulations and 
#' revised tests for dissociations. Neuropsychology, 19(3), 318.
#' 
#' @param caseData - A 1-by-2 matrix or data frame storing the data of the single 
#'                   case subject 
#' @param controlData - A n-by-2 matrix or data frame storing the data of 
#'                      control subjects for comparisions (n is the number of 
#'                      controls)
#' @param alphaLevel the significance level (alpha); default is 0.05
#' @return result - a data frame storing the t-value, degrees of freedom and p-value
#' 
#' Writtien by C.-H. Eric Chang, April, 2020
#' 
CrawfordRSDT <- function(caseData, controlData, alphaLevel = 0.05) {
  # Check input arguments
  # (1) whether caseData is a 1-by-2 matrix or data frame
  if (dim(caseData)[1]!=1) {
    stop(paste("Only 1 single case data is allowed. ", 
         "Number of rows for the single case data matrix/data frame should be 1.",
         sep = ""))
  } else if (dim(caseData)[2]!=2) {
    stop(paste("This function only compares performance of 2 tasks for a single case. ", 
               "Number of columns for the single case data matrix/data frame should be 2.",
               sep = ""))
  }
  # (2) whether controls data is a n-by-2 matrix or data frame
  if (dim(controlData)[2]!=2) {
    stop(paste("This function only compares performance of 2 tasks for a single case. ", 
               "Number of columns for the control subjects data matrix/data frame should be 2.",
               sep = ""))
  }
  # (3) whether there are missing valeus in the single case data
  if (sum(is.na(caseData))!=0) {
    stop("There are missing values in the single case data.")
  }
  # (4) whether there are issing valeus in the control subjects data
  if (sum(is.na(controlData))!=0) {
    cat("There are missing values in the control subjects data. Remove controls with missing values.","\n")
    controlData_ori <- controlData
    indNA <- which(rowSums(is.na(controlData))!=0) # which control(s) has missing values
    controlData <- controlData[-indNA,] # remove controls with missing values
  }
  # (5) if caseData and/or controlData is data frame, convert them to matrix
  if (is.data.frame(caseData)) {
    caseData_df <- caseData
    caseData <- data.matrix(caseData)
  }
  if (is.data.frame(controlData)) {
    controlData_df <- controlData
    controlData <- data.matrix(controlData)
  }
  
  # Parameters
  averC1 <- mean(controlData[,1]) # mean of task 1 in controls
  averC2 <- mean(controlData[,2]) # mean of task 2 in controls
  stdC1 <- sd(controlData[,1]) # standard deviation of task 1 in controls
  stdC2 <- sd(controlData[,2]) # standard deviation of task 2 in controls
  nC <- dim(controlData)[1] # number of control subjects
  r <- cor(controlData[,1],controlData[,2]) # correlation between task 1 and task 2 in controls
  y <- qt(1-alphaLevel/2, df = nC-1) # critical t value (2-tailed) on n-1 degrees of freedom
  a <- (1+r)*(1-r^2)
  b <- (1-r)*( (4*((nC-1)^2)) + (4*(1+r)*(nC-1)) + ((1+r)*(5+r)) )
  c <- -2*(( (caseData[1]-averC1)/stdC1 - (caseData[2]-averC2)/stdC2 )^2) * ((nC*((nC-1)^2)) / (nC+1))
  
  # Compute psi (formula 6 in the paper)
  psi <- ( (caseData[1]-averC1)/stdC1 - (caseData[2]-averC2)/stdC2 ) /
    sqrt( ((nC+1)/nC) * ( (2-2*r) +  ((2*(1-r^2))/(nC-1)) + (((5+y^2)*(1-r^2))/(2*((nC-1)^2))) + ((r*(1+y^2)*(1-r^2))/(2*((nC-1)^2))) ) ) 
  
  # Statistical test results
  h <- ifelse(abs(psi) > y, 1, 0) # whether null hypothesis is rejected (1:yes; 0:no)
  tVal <- sqrt( (-b + sqrt(b^2-4*a*c)) / (2*a) )
  tVal <- ifelse(psi>0, tVal, -tVal)
  pval <- 2*(1-pt(abs(tVal), df=nC-1)) #two-tailed p-value
  
  # Effect size (as a z-score)
  diff1 <- caseData[1]-averC1
  diff2 <- caseData[2]-averC2
  corVal <- cor(controlData[,1], controlData[,2]) # correlation between the 2 tasks in control
  effSize = ( (diff1/stdC1) - (diff2/stdC2) ) / sqrt(2-2*corVal)
  
  # Return the results
  result <- data.frame(t = tVal, df = nC-1, p = pval, Zdcc=effSize)
  return(result) 

  }