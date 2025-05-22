## load data objects to package
setwd("C:/Users/kdmulligan/OneDrive - University of Iowa/Dissertation/ABM/Data")
# getwd()
# list.files()
CDI_los_dist_mdc_tran6days <-
  readRDS("CDI_los_dist_mdc_tran6days.RDS")

noCDI_los_dist_mdc_tran6days <-
  readRDS("noCDI_los_dist_mdc_tran6days.RDS")


## add data to r package
use_data(CDI_los_dist_mdc_tran6days)
use_data(noCDI_los_dist_mdc_tran6days)
