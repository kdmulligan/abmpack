utils::globalVariables(c("susceptible", "latent", "asymptomatic", "incubation",
                         "symptomatic", "clin_res", "death"))
#### function 8: pat_in_one_state()

#' @title Verify patient(s) in one disease state
#'
#' @description Take a vector of patients IDs and check if patients are in
#' more than one state
#'
#' @param pats vector of patient ids of interest
#'
#' @return logical vector
#' @importFrom tibble tibble
#' @importFrom dplyr mutate if_else
#' @export

pat_in_one_state = function(pats) {
  states = tibble(
    'sus' = pats %in% (susceptible@i + 1),
    'lat' = pats %in% (latent@i + 1),
    'asy' = pats %in% (asymptomatic@i + 1),
    'inc' = pats %in% (incubation@i + 1),
    'sym' = pats %in% (symptomatic@i + 1),
    'clr' = pats %in% (clin_res@i + 1),
    'dea' = pats %in% (death@i + 1)
  )

  check_df = states |>
    mutate(
      sum = sus + lat + asy + inc + sym +clr + dea,
      check = if_else(sum == 1, TRUE, FALSE)
    )
  return(check_df$check)
}
