utils::globalVariables(c("hcup_id_vec"))

#' @title Update list of available rooms based on occupancy
#'
#' @description This function has no parameters because it uses the room_list
#' object from the environment. Additionally, it has no return object because
#' it replaces icu_avail and non_avail in the environment.
#' @param rm_list list object with room info
#'
#' @return Data frame with number of rows equal to the number of patient moves per day
#' @export

update_available_rooms = function(rm_list) {
  avail_idx = rm_list$rid_uniq[!(rm_list$rid_uniq) %in% (rm_list$occup@i + 1)]
  icu_avail_temp <- vector("list", length(hcup_id_vec))
  # non_avail_temp <- vector("list", length(hcup_id_vec))
  names(icu_avail_temp) <- as.character(hcup_id_vec)
  # names(non_avail_temp) <- as.character(hcup_id_vec)
  non_avail_temp <- icu_avail_temp
  for(j in 1:length(hcup_id_vec)) {
    # index: facility, ICU/NON, available
    icu_avail_temp[[j]] <- rm_list$rid_uniq[(rm_list$hcup_id == hcup_id_vec[j]) & (rm_list$icu == 1) & (rm_list$rid_uniq %in% avail_idx)]
    non_avail_temp[[j]] <- rm_list$rid_uniq[(rm_list$hcup_id == hcup_id_vec[j]) & (rm_list$icu == 0) & (rm_list$rid_uniq %in% avail_idx)]
  }
  ## original
  # assign("icu_avail", icu_avail_temp, envir = .GlobalEnv)
  # assign("non_avail", non_avail_temp, envir = .GlobalEnv)
  ## parent envir
  assign("icu_avail", icu_avail_temp, envir = parent.frame())
  assign("non_avail", non_avail_temp, envir = parent.frame())
}
