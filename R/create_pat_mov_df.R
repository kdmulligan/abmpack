#### function 5: create_pat_mov_df()

#' @title Construct data frame of patient moves for the day
#'
#' @description Take all the patients in the hospital (their id, visit LOS, &
#' current room type icu/non) and construct a data frame of their movements for
#' the coming day. The input vectors should all be of the same length, with length
#' equal to the number of patients in the hospital on that day. If a patient's
#' drawn_moves for the day are "0,0" they will not be in the resulting data set.
#'
#'
#' @param ids vector of patient ids
#' @param start_loc vector of patient starting room location (1/0) for (icu/non) for the day
#' @param drawn_moves vector of the number of icu/non moves for the day, designed to take
#'  output from my `sample_day_mvts()` rcpp function
#'
#' @return Data frame with number of rows equal to the number of patient moves per day
#'
#' @importFrom dplyr case_when
#' @export

# sourceCpp("rcpp/get_day_mvts.cpp")

create_pat_mov_df <- function(ids, start_loc, drawn_moves) {
  idx_keep = which(drawn_moves != "0,0")
  ids_keep = ids[idx_keep]
  icu_n = as.numeric(substr(drawn_moves[idx_keep], 1, 1))
  non_n = as.numeric(substr(drawn_moves[idx_keep], 3, 3))
  end_loc_keep = case_when(
    # prob end in ICU if start in ICU
    drawn_moves[idx_keep] %in% c("1,1", "1,2", "2,1") & start_loc[idx_keep] == 1 ~ rbinom(1, 1, 0.511),
    # prob end in ICU if start in NON
    drawn_moves[idx_keep] %in% c("1,1", "1,2", "2,1") & start_loc[idx_keep] == 0 ~ rbinom(1, 1, 1 - 0.92),
    # no randomness in end location
    drawn_moves[idx_keep] %in% c("1,0", "2,0", "3,0") ~ 1,
    drawn_moves[idx_keep] %in% c("0,1", "0,2", "0,3") ~ 0,
  )
  ## pass to rcpp function
  df <- get_day_mvts_cpp(
    ids = ids_keep, icu_num = icu_n,
    non_num = non_n, end_loc = end_loc_keep
  )
  return(df)
}
