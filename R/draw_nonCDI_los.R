utils::globalVariables(c("noCDI_los_dist_mdc_tran6days"))
#### function 3: draw_nonCDI_los()

#' @title draw non-CDI los before simulation
#'
#' @description If a patient has CDI according to the HCUP data, their LOS
#' needs to be redrawn from the non_CDI LOS. The new LOS is drawn from the
#' matching transfer/MDC distribution with days between 1-365
#'
#' @param dist a list object with los distributions by mdc and transfer type
#' @param transfer transfer type of the hospitalization: not a transfer,
#' transfer not last, transfer last
#' @param md_cat major diagnostic category, 0-25
#'
#' @return Returns a integer (or vector, function is vectorized) with the new LOS
#'
#' @importFrom dplyr case_when
#' @export

draw_nonCDI_los <- function(dist = noCDI_los_dist_mdc_tran6days, transfer, md_cat){
  ## distribution: noCDI_los_dist_mdc_tran
  transfer = case_when(
    transfer == 1 ~ "non_tran",
    transfer == 2 ~ "tran_last",
    transfer == 3 ~ "tran_not_last"
  )
  out <- mapply(function(arg1, arg2)
    sample(
      x = 1:365,
      size = 1,
      prob = dist[[arg1]][[paste0("mdc_", arg2)]]),
    arg1 = transfer,
    arg2 = md_cat
  )
  return(out)
}

