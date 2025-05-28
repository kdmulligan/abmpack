utils::globalVariables(
  c("age", "key", "key1", "key2", "key3", "key4", "age.x", "age.y",
    "pointoforiginub04", "pointoforiginub04.x", "pointoforiginub04.y",
    "adrgsev.x", "adrgsev.y", "age_cat", "adrgriskmortality.x", "adrgriskmortality.y"))
#### function 1: add_supp_hcup_dat()
#' @title combine additional HCUP data with data originally pulled from HCUP
#'
#' @description The original HCUP data took a while to run and later I decided
#' there was additional data I wanted to pull from HCUP so it was faster to
#' pull it separately and then join it based on the key. This function takes
#' the original data (`base_dat`) and joins the supplementary data
#' to it (`supp_dat`).
#'
#' @param base_dat original processed data from hcup
#' @param supp_dat additional/supplementary processed data from HCUP
#'
#' @return Returns a tibble with the original data and the new data combined.
#'
#' @importFrom dplyr rename filter left_join mutate relocate select rowwise ungroup bind_rows
#' @importFrom tidyr separate_wider_delim replace_na
#' @importFrom stringr str_detect str_c
#' @export

add_supp_hcup_dat <- function(base_dat, supp_dat) {
  ## data with internal transfers (>1 key)
  dat_inter_tran <-
    base_dat |>
    rename(age_cat = age) |>
    filter(str_detect(key, " "))
  ## data with no internal transfers (only 1 key)
  dat_no_inter_tran <-
    base_dat |>
    rename(age_cat = age) |>
    filter(!str_detect(key, " "))
  ## add to transfers data first
  tran <-
    dat_inter_tran |>
    separate_wider_delim(
      cols = key, delim = " ",
      names = c("key1", "key2", "key3", "key4"),
      too_few = "align_start",
      too_many = "drop" ## 105 segments have > 4 transfers, UT
    ) |>
    left_join(supp_dat, join_by(key1 == key)) |>
    left_join(supp_dat, join_by(key2 == key)) |>
    mutate(
      age = if_else(!is.na(age.x), age.x, age.y),
      pointoforiginub04 = if_else(!is.na(pointoforiginub04.x), pointoforiginub04.x, pointoforiginub04.y),
      tran_fr_snf = if_else(pointoforiginub04 == "5" & age > 0, 1, 0),
      adrgriskmortality = if_else(!is.na(adrgriskmortality.x), adrgriskmortality.x, adrgriskmortality.y),
      adrgsev = if_else(!is.na(adrgsev.x), adrgsev.x, adrgsev.y),
      .before = age.x
    ) |>
    relocate(age, .after = age_cat) |>
    dplyr::select(-c(age.x, age.y, pointoforiginub04.x, pointoforiginub04.y, adrgriskmortality.x, adrgriskmortality.y, adrgsev.x, adrgsev.y)) |>
    replace_na(list(key3 = "", key4 = "")) |>
    rowwise() |>
    mutate(
      key = str_c(key1, key2, key3, key4, sep = " "), .before = key1
    ) |>
    ungroup() |>
    dplyr::select(-c(key1, key2, key3, key4))

  ## add to non-transfer data
  notran <-
    dat_no_inter_tran |>
    left_join(supp_dat, join_by(key)) |>
    relocate(age, .after = age_cat) |>
    mutate(
      tran_fr_snf = if_else(pointoforiginub04 == "5" & age > 0, 1, 0),
      .after = pointoforiginub04
    )
  ## combine datasets
  dat <- bind_rows(notran, tran)
  return(dat)

}
