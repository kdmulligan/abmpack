utils::globalVariables(c("hcupids_vec"))

#' @title Update list of available rooms based on occupancy
#'
#' @description This function has no parameters because it uses the room_list
#' object from the environment. Additionally, it has no return object because
#' it replaces icu_avail and non_avail in the environment.
#' @param rm_list list object with room info
#' @param hcupids_vec vector of hcup ids
#'
#' @return Data frame with number of rows equal to the number of patient moves per day
#' @export

update_available_rooms = function(rm_list, hcupids_vec) {
  avail_idx = rm_list$rid_uniq[!(rm_list$rid_uniq) %in% (rm_list$occup@i + 1)]
  icu_avail_temp <- vector("list", length(hcupids_vec))
  # non_avail_temp <- vector("list", length(hcupids_vec))
  names(icu_avail_temp) <- as.character(hcupids_vec)
  # names(non_avail_temp) <- as.character(hcupids_vec)
  non_avail_temp <- icu_avail_temp
  for(j in 1:length(hcupids_vec)) {
    # index: facility, ICU/NON, available
    icu_avail_temp[[j]] <- rm_list$rid_uniq[(rm_list$hcup_id == hcupids_vec[j]) & (rm_list$icu == 1) & (rm_list$rid_uniq %in% avail_idx)]
    non_avail_temp[[j]] <- rm_list$rid_uniq[(rm_list$hcup_id == hcupids_vec[j]) & (rm_list$icu == 0) & (rm_list$rid_uniq %in% avail_idx)]
  }
  ## original
  assign("icu_avail", icu_avail_temp, envir = .GlobalEnv)
  assign("non_avail", non_avail_temp, envir = .GlobalEnv)
  ## parent envir
  # assign("icu_avail", icu_avail_temp, envir = parent.frame())
  # assign("non_avail", non_avail_temp, envir = parent.frame())
}
