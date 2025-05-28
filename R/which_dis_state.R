utils::globalVariables(c("susceptible", "latent", "asymptomatic", "incubation",
                         "symptomatic", "clin_res", "death"))
#### function 7: which_dis_state()

#' @title check disease states from a vector of patient IDs
#'
#' @description take vector of patient IDS and return disease status of each patient
#'
#' @param pats vector of patient ids of interest
#'
#' @return character vector of disease states

which_dis_state = function(pats) {
  which_sus = which(pats %in% (susceptible@i + 1))
  which_lat = which(pats %in% (latent@i + 1))
  which_asy = which(pats %in% (asymptomatic@i + 1))
  which_inc = which(pats %in% (incubation@i + 1))
  which_sym = which(pats %in% (symptomatic@i + 1))
  which_clr = which(pats %in% (clin_res@i + 1))
  which_dea = which(pats %in% (death@i + 1))

  elements = c(
    rep("suscep.", length(which_sus)), rep("latent", length(which_lat)),
    rep("asymp.", length(which_asy)), rep("incub.", length(which_inc)),
    rep("symp.", length(which_sym)), rep("clinres", length(which_clr)),
    rep("death", length(which_dea))
  )
  out = elements[order(c(which_sus, which_lat, which_asy, which_inc, which_sym, which_clr, which_dea))]
  return(out)
}
