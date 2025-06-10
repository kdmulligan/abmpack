utils::globalVariables(
  c("CDI_los_dist_mdc_tran6days", "hospid_max", "hcup_id", "max_los"))
#### function 2: draw_CDI_los()
#' @title draw CDI los based on current day during simulation
#'
#' @description If a patient gets CDI during the simulation, their LOS should
#' be adjusted to be longer. Draw a new LOS from the matching transfer/MDC
#' distribution truncated above the current LOS + 1. If the current_day is
#' great than or equal to the maximum los for the distribution, then the new LOS
#' is not drawn from the truncated distribution, instead it is choosen from the
#' maximum of the (current_day + 1) or los_hcup
#'
#' @param dist a list object with los distributions by mdc and transfer type
#' @param current_day day of hospitalization for the patient
#' @param transfer transfer type of the hospitalization: not a transfer,
#' transfer not last, transfer last
#' @param md_cat major diagnostic category, 0-25
#' @param hcup_los the patient's los in the HCUP data
#' @param hosp_id hospital the patient is from
#'
#' @return Returns a integer (or vector, function is vectorized) as the new LOS
#'
#' @importFrom dplyr filter pull case_when
#' @export
#'
# CDI_los_dist_mdc_tran6days <- readRDS(paste0(getwd(), "/Data/CDI_los_dist_mdc_tran6days.RDS"))

draw_CDI_los <- function(dist = CDI_los_dist_mdc_tran6days, current_day, transfer, md_cat, hcup_los, hosp_id) {
  ## distribution: CDI_los_dist_mdc_tran
  max_los_facil <- hospid_max |>
    filter(hcup_id == hosp_id) |>
    pull(max_los)
  if(current_day == 0 & transfer == 0) {
    transfer = 1
    md_cat = 6
  }
  transfer = case_when(
    transfer == 1 ~ "non_tran",
    transfer == 2 ~ "tran_last",
    transfer == 3 ~ "tran_not_last"
  )
  dist_max <- max(which(dist[[transfer]][[paste0("mdc_", md_cat)]] != 0))
  max_to_use <- min(max_los_facil, dist_max)

  out <- mapply(
    function(arg1, arg2, arg3, arg4, arg5)
      if ((arg1 + 1 == arg5) &  arg1 < (max(which(dist[[arg2]][[paste0("mdc_", arg3)]] != 0)) - 1)) {
        arg5
      } else if (arg1 < (max(which(dist[[arg2]][[paste0("mdc_", arg3)]] != 0)) - 1)) {
        print(paste0("current day: ", arg1, "; transfer: ", arg2, "; mdc: ", arg3,
                     "; hcup_los: ", arg4, "; min facility/dist: ", arg5
                     ))
        sample(
          x = (arg1 + 1):arg5,
          size = 1,
          prob = dist[[arg2]][[paste0("mdc_", arg3)]][(arg1 + 1):arg5]
        )
      } else {
        max(arg1 + 1, arg4)
      }
    ,
    arg1 = current_day,
    arg2 = transfer,
    arg3 = md_cat,
    arg4 = hcup_los,
    arg5 = max_to_use
  )
  return(out)
}
