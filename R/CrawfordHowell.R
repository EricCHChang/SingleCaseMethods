#' Crawford-Howell t-test for single case studies
#' Compare the score of the single case to controls
#' The formula is based on Crawford & Howell (1998) 
#' Comparing an individual's test score against norms derived from small samples. 
#' The Clinical Neuropsychologist, 12(4), 482-486.
#' 
#' @param case Data of a single case subject 
#' @param control Data of a group of subjects for comparison 
#' 
#' Writtien by C.-H. Eric Chang, Jan, 2019
#' 
CrawfordHowell <- function(case, control){
  tval <- (case - mean(control)) / (sd(control)*sqrt((length(control)+1) / length(control)))
  degfree <- length(control)-1
  pval <- 2*(1-pt(abs(tval), df=degfree)) #two-tailed p-value
  effSize <- (case - mean(control)) / sd(control)
  perct_ctrlBelowCase = pnorm(effSize, 0, 1)
  result <- data.frame(t = tval, df = degfree, p=pval, Zcc=effSize, perct_ctrlBelowCase=perct_ctrlBelowCase)
  return(result)
}