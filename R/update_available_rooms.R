utils::globalVariables(c("room_list", "hcup_id_vec"))
#### function 6: update_available_rooms()

#' @title Update list of available rooms based on occupancy
#'
#' @description This function has no parameters because it uses the room_list
#' object from the environment. Additionally, it has no return object because
#' it replaces icu_avail and non_avail in the environment.
#'
#' @return Data frame with number of rows equal to the number of patient moves per day

update_available_rooms = function() {
  avail_idx = room_list$rid_uniq[!(room_list$rid_uniq) %in% (room_list$occup@i + 1)]
  icu_avail_temp <- vector("list", length(hcup_id_vec))
  # non_avail_temp <- vector("list", length(hcup_id_vec))
  names(icu_avail_temp) <- as.character(hcup_id_vec)
  # names(non_avail_temp) <- as.character(hcup_id_vec)
  non_avail_temp <- icu_avail_temp
  for(j in 1:length(hcup_id_vec)) {
    # index: facility, ICU/NON, available
    icu_avail_temp[[j]] <- room_list$rid_uniq[(room_list$hcup_id == hcup_id_vec[j]) & (room_list$icu == 1) & (room_list$rid_uniq %in% avail_idx)]
    non_avail_temp[[j]] <- room_list$rid_uniq[(room_list$hcup_id == hcup_id_vec[j]) & (room_list$icu == 0) & (room_list$rid_uniq %in% avail_idx)]
  }
  ## original
  # assign("icu_avail", icu_avail_temp, envir = .GlobalEnv)
  # assign("non_avail", non_avail_temp, envir = .GlobalEnv)
  ## parent envir
  assign("icu_avail", icu_avail_temp, envir = parent.frame())
  assign("non_avail", non_avail_temp, envir = parent.frame())
}
