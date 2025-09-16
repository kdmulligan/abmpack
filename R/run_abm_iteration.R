utils::globalVariables(
  c("hospid_max", "rooms", "uniq_pat_dat", "hcws", "hcup", "aday_sim", "hcup_id",
    "adrgriskmortality", "patid", "los_rem_frday0", "days2next", "class_sim", "viz_num",
    "subseq_tran_seg", "viz_key", "n", "risk_1", "risk_2", "risk_3", "risk_4", "fid", "time_block",
    "room_type", "hcw_mvts", "day_counter", "rid", "hid", "total_sec", "rid_uniq", "hid_ss",
    "hid_uniq", 'prop_soap_in', "prop_soap_out", "prob_contam_col", "prob_contam_contam",
    "prob_contam_fr_patdis", "prob_room_nocontam_hcw_col", "prob_room_nocontam_hcw_contam",
    "prob_contam", "prob_hcw_now_contam", "room_list", "hcup_id_vec"
  )
)

#' @title Run one ABM simulation
#' @description Take a input of parameters and indicators to run simulation by
#'
#' @param n_days number of days for the simulation
#' @param t_symp_room transmission from symp pat to room prior value
#' @param t_hcw_room transmission from hcw to room prior value
#' @param t_room_hcw transmission from room to hcw prior value
#' @param t_room_patge65 transmission from room to patient 65 or older prior value
#' @param t_prob_room_contam_init probability a room is contam prior value
#' @param l_binom_latent binomial distribution probability for latent period prior value
#' @param ind_col_pat_trans indicator for colonized patient transmission
#' @param ind_asymp_pat_trans indicator for asymptomatic patient transmission
#' @param ind_asymp_hcw_trans indicator for asymptomatic hcw transmission
#' @param ind_recur_trans indicator for recurrent patient transmission
#' @param ind_transfer_pat_trans indicator for transfer patient transmission
#' @param ind_revisit_pat_trans indicator for revisit patient tranmission
#' @param revisit_n_days number of days to consider for a revisit. Must be less
#' than or equal to 6. Only important to set if setting `ind_revisit_pat_trans` to zero.
#' @param ind_symp_trans indicator for symptomatic patient transmission
#' @param SEED seed for the simulation
#'
#' @return numeric value with total number of observed cases
#' @importFrom Matrix Matrix
#' @importFrom tibble tibble
#' @importFrom utils str
#' @importFrom dplyr mutate if_else filter arrange desc bind_rows left_join pull
#' select join_by inner_join ungroup group_by count slice_sample rename summarise
#' summarize
#' @importFrom tidyr replace_na pivot_wider
#' @importFrom stats rbinom runif rweibull rpois
#' @importFrom truncnorm rtruncnorm
#' @importFrom triangle rtriangle
#' @export

run_abm_iteration <- function(n_days = 72,
                              t_symp_room = 0.0005, t_hcw_room = 0.0005,
                              t_room_hcw = 0.0005, t_room_patge65 = 0.0005,
                              t_prob_room_contam_init = 0.005, l_binom_latent = 1,
                              ind_col_pat_trans = 1, ind_asymp_pat_trans = 1,
                              ind_asymp_hcw_trans = 1, ind_recur_trans = 1,
                              ind_transfer_pat_trans = 1, ind_revisit_pat_trans = 1,
                              revisit_n_days = 6, ind_symp_trans = 1, SEED = 1212
                              )
{
  set.seed(SEED)
  ## number of sim days
  total_days_sim = n_days
  ## parameter values
  tau_symp_room = t_symp_room
  tau_hcw_room = t_hcw_room
  tau_room_hcw = t_room_hcw
  tau_room_patge65 = t_room_patge65
  tau_prob_room_contam_init = t_prob_room_contam_init
  ## research question indicators
  ## RQ 1
  incl_col_pat_trans = ind_col_pat_trans
  incl_asymp_pat_trans = ind_asymp_pat_trans
  incl_asymp_hcw_trans = ind_asymp_hcw_trans
  ## RQ 2
  incl_recur_trans = ind_recur_trans
  ## RQ 3
  incl_transfer_pat_trans = ind_transfer_pat_trans
  ## RQ 4
  incl_revisit_pat_trans = ind_revisit_pat_trans
  rq_tran_days_bn = revisit_n_days          ## transfer def is interviz window of <= 6 days by default
  ## RQ 5
  incl_symp_trans = ind_symp_trans

  ## create list to hold results
  results = list(
    tot_symp_fr_incu = 0,
    tot_symp_bc_recur = 0,
    tot_death = 0,
    tot_clres = 0,
    tot_asymp = 0
  )
  obs_results = list(
    obs_tot_symp = 0,
    obs_recur = 0,
    obs_hcw_asymp = 0,
    obs_all_symp = 0,
    prev_feb = 0,
    prev_mar = 0,
    prop_rm_contam_end = 0
  )
  tot_pat_newsymp_day = vector(mode = "numeric", length = n_days)

  ## step 0: INITIALIZATION ####################################################

  ## tran_seg note: 1 = non_tran, 2 = tran_last, 3 = tran_not_last
  ## transmission ratio by age group compared to ge65
  age_lambdas = c(0.969, 0.286, 0.406, 0.619, 0.801, 1)
  # max los by hcup facility
  # hospid_max <-
  #   read_csv("Data/hcup_final_sim_2010dat.csv", show_col_types = FALSE) |>
  #   group_by(hcup_id) |>
  #   summarize(max_los = max(los_sim, na.rm = TRUE))
  ## facility ids HCUP
  hcup_id_vec <- unique(hospid_max$hcup_id)
  ## rooms list
  n_rooms = nrow(rooms)
  room_list = list(
    # static
    hcup_id = rooms$hcup_id,
    rid_uniq = rooms$rid_uniq,
    rid_ss = rooms$rid,
    icu = rooms$icu_flag,
    # dynamic
    contam = Matrix(0L, n_rooms, 1),
    contam_next = Matrix(0L, n_rooms, 1),
    occup = Matrix(0L, n_rooms, 1)
  )
  # print(str(room_list))
  ## patient list
  n_pat = nrow(uniq_pat_dat)
  pat_list = list(
    # static
    age_cat = uniq_pat_dat$age_cat_n,
    ge65_flag = if_else(uniq_pat_dat$age_cat_n == 6, 1, 0),
    # dynamic
    los_sim = Matrix(0L, n_pat, 1),
    days_rem = Matrix(0L, n_pat, 1),
    tran_stat = Matrix(0L, n_pat, 1),
    tran_sub_fac_diff = Matrix(0L, n_pat, 1),
    tran_days_bn = Matrix(0L, n_pat, 1),
    mdc = Matrix(0L, n_pat, 1),
    facility = Matrix(0L, n_pat, 1),
    risk_mor = Matrix(0L, n_pat, 1),
    viz_key = Matrix(0L, n_pat, 1),
    days_to_next_viz = Matrix(0L, n_pat, 1),
    viz_num = Matrix(0L, n_pat, 1),
    new_symp_viz = Matrix(0L, n_pat, 1)
  )
  ge65_uniq_pat_idx = which(pat_list$ge65_flag == 1)

  ## hcw list
  n_hcw = nrow(hcws)
  hcw_list = list(
    # static
    hcup_id = hcws$hcup_id,
    job_cat = hcws$job_cat,
    prob_wash_in = hcws$prop_soap_in,
    prob_wash_out = hcws$prop_soap_out,
    # dynamic
    colonized = Matrix(0L, n_hcw, 1),
    contam = Matrix(0L, n_hcw, 1),
    contam_next = Matrix(0L, n_hcw, 1)
  )

  ## create sparse Matrices to hold compartmental disease model states
  susceptible = Matrix(0L, n_pat, 1)
  latent = Matrix(0L, n_pat, 1)
  asymptomatic = Matrix(0L, n_pat, 1)
  incubation = Matrix(0L, n_pat, 1)
  symptomatic = Matrix(0L, n_pat, 1)
  clin_res = Matrix(0L, n_pat, 1)
  after_clin_res_location = Matrix(0L, n_pat, 1)
  death = Matrix(0L, n_pat, 1)

  ## vectors for recurrent indices
  recur_1_idx = numeric(0)  ## indices of patients with 1 recurrence
  recur_2p_idx = numeric(0) ## indices of patients with 2+ recurrence
  all_admit_keys = numeric(0)

  ## data frame for queued patients
  queue_viz_keys_df = tibble(
    patid = numeric(),
    viz_key = numeric(),
    hcup_id = numeric(),
    adrgriskmortality = numeric(),
    los_sim = numeric()
  )  ## vector of viz keys for patients in queue

  ## lists of rooms per facility (all available at initialization here)
  icu_avail <- vector("list", length(hcup_id_vec))
  non_avail <- vector("list", length(hcup_id_vec))
  names(icu_avail) <- as.character(hcup_id_vec)
  names(non_avail) <- as.character(hcup_id_vec)
  for (j in 1:length(hcup_id_vec)) {
    # index: facility, ICU/NON
    icu_avail[[j]] <- room_list$rid_uniq[(room_list$hcup_id == hcup_id_vec[j]) &
                                           (room_list$icu == 1)]
    non_avail[[j]] <- room_list$rid_uniq[(room_list$hcup_id == hcup_id_vec[j]) &
                                           (room_list$icu == 0)]
  }

  ## set probability of hcws to be colonized upon initialization
  hcw_colon = rbinom(
    n_hcw,
    size = 1,
    prob = rtriangle(
      n = n_hcw,
      a = 0.001,
      b = 0.05,
      c = 0.025
    )
  )
  hcw_colon = which(hcw_colon == 1)
  hcw_list$colonized[hcw_colon] = 1L
  ## set probability of hcws to be contam upon initialization
  # tau_prob_room_contam_init = 0.005
  probs = 1 - (0.16^(hcws$avg_rms_day * tau_prob_room_contam_init))
  hcw_contam = rbinom(n_hcw, size = 1, prob = probs)
  hcw_contam = which(hcw_contam == 1)
  # length(hcw_contam) / n_hcw
  hcw_list$contam[hcw_contam] = 1L
  ## set probability of rooms contam upon initialization
  # tau_prob_room_contam_init
  room_contam = rbinom(n = n_rooms, size = 1, prob = tau_prob_room_contam_init)
  room_contam = which(room_contam == 1)
  room_list$contam[room_contam] = 1L

  ## admit patients in hosp on day 0
  to_admit_df <- hcup |>
    filter(aday_sim == 0)

  # assign rooms
  new_rooms <- vector("list", length(hcup_id_vec))
  for (i in 1:length(hcup_id_vec)) {
    ## this for loop is fast
    new_rooms[[i]] <-
      to_admit_df |>
      filter(hcup_id == hcup_id_vec[i]) |>
      arrange(desc(adrgriskmortality)) |>
      dplyr::select(patid, adrgriskmortality) |>
      rcpp_assign_rooms_cpp_seed(
        pat_risks = _,
        icu = icu_avail[[as.character(hcup_id_vec[i])]],
        non = non_avail[[as.character(hcup_id_vec[i])]],
        seed = SEED
      )
  }
  new_rooms <- bind_rows(new_rooms)
  # set rooms as occupied
  room_list$occup[new_rooms$assigned_room] = new_rooms$patid
  # add patient inform?ation for new hcup admits
  all_admit_keys = c(all_admit_keys, to_admit_df$viz_key) ## dynamic vector of all used keys
  # admit patients
  to_admit_pat = to_admit_df$patid
  pat_list$los_sim[to_admit_pat] = to_admit_df$los_sim
  pat_list$days_rem[to_admit_pat] =
    to_admit_df |>
    dplyr::select(patid) |>
    left_join(uniq_pat_dat |> select(patid, los_rem_frday0),
              join_by(patid)) |>
    pull(los_rem_frday0)
  pat_list$tran_stat[to_admit_pat] = to_admit_df$tran_seg
  pat_list$tran_sub_fac_diff[to_admit_pat] = to_admit_df$fac_diff_flag
  pat_list$tran_days_bn[to_admit_pat] = to_admit_df$days2next
  pat_list$mdc[to_admit_pat] = to_admit_df$mdc
  pat_list$facility[to_admit_pat] = to_admit_df$hcup_id
  pat_list$risk_mor[to_admit_pat] = to_admit_df$adrgriskmortality
  pat_list$viz_key[to_admit_pat] = to_admit_df$viz_key
  pat_list$viz_num[to_admit_pat] = to_admit_df$viz_num

  next_viz_dat = to_admit_df |>
    filter(!is.na(days2next))
  pat_list$days_to_next_viz[next_viz_dat$patid] = next_viz_dat$days2next

  # disease status for initial pat, find pat with CDI visit
  cdi_admits = to_admit_df |>
    select(patid) |>
    left_join(uniq_pat_dat |> select(patid, class, class_sim),
              join_by(patid)) |>
    filter(class_sim == "CDI")
  ## patients entering at symptomatic
  idx_to_symp = cdi_admits |>
    filter(class %in% c("R", "HAI", "T")) |>
    pull(patid)
  symptomatic[idx_to_symp] =
    ceiling(triangle::rtriangle(
      n = length(idx_to_symp),
      a = 3,
      b = 12,
      c = 8
    ) *
      runif(
        n = length(idx_to_symp),
        min = 0.1,
        max = 0.9
      )) ## runif gets proportion of state days that occured before day 0

  ## patients entering at incubation
  idx_to_incub = cdi_admits |>
    filter(class == "CAI") |>
    pull(patid)
  incubation[idx_to_incub] =
    ceiling(triangle::rtriangle(
      n = length(idx_to_incub),
      a = 1,
      b = 10,
      c = 4
    ) * runif(
      n = length(idx_to_incub),
      min = 0.1,
      max = 0.9
    ))
  ## non CDI patients
  ## set probability of noCDI patients to be asymp at first viz
  maybe_asymp_pat = to_admit_df$patid[!to_admit_df$patid %in% cdi_admits$patid]
  idx_pat_asymp = maybe_asymp_pat[which(rbinom(length(maybe_asymp_pat), size = 1, prob = 0.00785) == 1)]
  asymptomatic[idx_pat_asymp] = ceiling(rweibull(length(idx_pat_asymp), 1, 43) * runif(
    n = length(idx_pat_asymp),
    min = 0.1,
    max = 0.9
  ))
  idx_pat_sus = maybe_asymp_pat[!maybe_asymp_pat %in% idx_pat_asymp]
  susceptible[idx_pat_sus] = 1L

  ## first patients that will be discharged
  idx_to_discharge = uniq_pat_dat |> filter(los_rem_frday0 == 0) |> pull(patid)

  ## set current obs symp cdi
  symp_pat_idx = symptomatic@i + 1
  symp_in_hosp = symp_pat_idx[symp_pat_idx %in% (room_list$occup@x)]
  symp_in_hosp_yesterday = symp_in_hosp
  #
  recur_in_hosp_yesterday = c(recur_1_idx, recur_2p_idx)

  ### start simulation after initalizations
  start_time = Sys.time()
  for (d in 1:total_days_sim) {
      print(paste0("start of DAY ", d))
      # print(c(452, 520, 1440, 2498) %in% room_list$occup@x)
      #       print(c(452, 520, 1440, 2498) %in% idx_to_discharge)
      # print(paste0("pat 46184 dis stat: ", which_dis_state(46184)))
      # print(paste0("patient 55718 state:", which_dis_state(55718)))
      # print(paste0("patient 55718 days rem: ", pat_list$days_rem[55718]))
      # print(paste0("patient 55718 in the hospital? ", 55718 %in% room_list$occup@x))


      ## step 1: pat days in disease state ###########################################
      ## update patient days since current disease state
        latent_end = (latent@i + 1)[which(latent@x == 1)]
        latent[latent_end] = 0L
        latent@x = latent@x - 1

        asymptomatic_end = (asymptomatic@i + 1)[which(asymptomatic@x == 1)]
        asymptomatic[asymptomatic_end] = 0L
        asymptomatic@x = asymptomatic@x - 1

        incubation_end = (incubation@i + 1)[which(incubation@x == 1)]
        incubation[incubation_end] = 0L
        incubation@x = incubation@x - 1

        symptomatic_end = (symptomatic@i + 1)[which(symptomatic@x == 1)]
        symptomatic[symptomatic_end] = 0L
        symptomatic@x = symptomatic@x - 1

        clin_res_end = (clin_res@i + 1)[which(clin_res@x == 1)]
        clin_res[clin_res_end] = 0L
        clin_res@x = clin_res@x - 1
      ## step 2: pat disease comparts ################################################
      ## update patient disease compartments if applicable

      n_to_move = length(latent_end) + length(asymptomatic_end) +
        length(incubation_end) + length(symptomatic_end) + length(clin_res_end)

      if (n_to_move > 0) {
        ## tran_1: from latent period to asymptotic infection indices
        # time in asymptotic infection: weibull(shape = 1, scale = 43) + 1
        to_asymp_idx = latent_end[which(rbinom(
          n = length(latent_end),
          size = 1,
          prob = 0.9
        ) == 1)]
        if (length(to_asymp_idx) > 0) {
          asymptomatic[to_asymp_idx] = round(rtruncnorm(
            n = length(to_asymp_idx),
            a = 14,
            b = 100,
            mean = 41.5,
            sd = 22.1
          ),
          0)
          # asymptomatic[to_asymp_idx] = round(rweibull(length(to_asymp_idx), 1, 43), 0) + 1
          results$tot_asymp = results$tot_asymp + length(to_asymp_idx)
        }
        ## tran_2: from latent period to incubation period
        # time in incubation period: triangle(6, 3, 9)
        to_incub_idx = latent_end[!latent_end %in% to_asymp_idx]
        if (length(to_incub_idx) > 0) {
          incubation[to_incub_idx] = round(triangle::rtriangle(
            n = length(to_incub_idx),
            a = 1,
            b = 10,
            c = 4
          ))
        }
        ## tran_3: from incubation period to symptomatic infection (everyone)
        ## equals draws incubation period number of days
        if (length(incubation_end) > 0) {
          symptomatic[incubation_end] = round(triangle::rtriangle(
            n = length(incubation_end),
            a = 5,
            b = 15,
            c = 10
          ), 0)
          results$tot_symp_fr_incu = results$tot_symp_fr_incu + length(incubation_end)
        }
        to_death = numeric()
        if (length(symptomatic_end) > 0) {
          ## tran_4: from symptomatic infection to death
          eval_logic = symptomatic_end %in% ge65_uniq_pat_idx
          symptomatic_end_ge65_idx = symptomatic_end[eval_logic]
          symptomatic_end_l65_idx = symptomatic_end[!eval_logic]
          to_death = c(symptomatic_end_ge65_idx[which(rbinom(
            n = length(symptomatic_end_ge65_idx),
            size = 1,
            prob = 0.091
          ) == 1)],
          symptomatic_end_l65_idx[which(rbinom(
            n = length(symptomatic_end_l65_idx),
            size = 1,
            prob = 0.0244
          ) == 1)])
          if (length(to_death) > 0) {
            death[to_death] = 1L
            results$death = results$death + length(to_death)
          }
          ## tran_5: from symptomatic infection to clincially resolved
          ## patients that do not die move to clinically resolved
          to_clin_res = symptomatic_end[!symptomatic_end %in% to_death]
          results$tot_clres = results$tot_clres + length(to_clin_res)
          ### patients less than 65 to susceptible or recurrent based on # of previous recur.
          to_clin_res_ge65_idx = to_clin_res[to_clin_res %in% ge65_uniq_pat_idx]
          ## >=65 with no prior recurrences
          no_recur_idx = to_clin_res_ge65_idx[!(
            to_clin_res_ge65_idx %in% recur_2p_idx |
              to_clin_res_ge65_idx %in% recur_1_idx
          )]
          one_recur_idx = to_clin_res_ge65_idx[to_clin_res_ge65_idx %in% recur_1_idx]
          two_plus_recur_ids = to_clin_res_ge65_idx[to_clin_res_ge65_idx %in% recur_2p_idx]
          ## indices for patients going to REC after CR
          recur_after_cr_ge65 = c(no_recur_idx[which(rbinom(
            n = length(no_recur_idx),
            size = 1,
            prob = 0.347
          ) == 1)], one_recur_idx[which(rbinom(
            n = length(one_recur_idx),
            size = 1,
            prob = 0.597
          ) == 1)], two_plus_recur_ids[which(rbinom(
            n = length(two_plus_recur_ids),
            size = 1,
            prob = 0.584
          ) == 1)])
          ## indices for patients going to SUSceptible after CR
          sus_after_cr_ge65 = to_clin_res_ge65_idx[!to_clin_res_ge65_idx %in% recur_after_cr_ge65]
          ### patients less than 65 to susceptible or recurrent
          to_clin_res_l65_idx = to_clin_res[!to_clin_res %in% ge65_uniq_pat_idx]
          recur_after_cr_l65 = to_clin_res_l65_idx[which(rbinom(
            n = length(to_clin_res_l65_idx),
            size = 1,
            prob = 0.1311
          ) == 1)]
          sus_after_cr_l65 = to_clin_res_l65_idx[!to_clin_res_l65_idx %in% recur_after_cr_l65]

          ## patients going to CR then susceptible
          n_cr_sus = (length(sus_after_cr_ge65) + length(sus_after_cr_l65))
          if (n_cr_sus > 0) {
            clin_res[c(sus_after_cr_ge65, sus_after_cr_l65)] =
              round(rtruncnorm(
                n = n_cr_sus,
                a = 7,
                b = 90,
                mean = 29.7,
                sd = 11.3
              ))
          }
          ## patients going to CR then recurrent infection (symptomtic)
          n_cr_rec = (length(recur_after_cr_ge65) + length(recur_after_cr_l65))
          if (n_cr_rec > 0) {
            clin_res[c(recur_after_cr_ge65, recur_after_cr_l65)] =
              round(triangle::rtriangle(
                n = n_cr_rec,
                a = 14,
                b = 96,
                c = 18
              ))
          }

          after_clin_res_location[c(sus_after_cr_ge65, sus_after_cr_l65)] = 7L
          after_clin_res_location[c(recur_after_cr_ge65, recur_after_cr_l65)] = 6L
        }

        if (length(clin_res_end) > 0) {
          ## tran_6: from clinically resolved to symptomatic infection (recurrent)
          indices_r = clin_res_end[clin_res_end %in% ((after_clin_res_location@i + 1)[which(after_clin_res_location@x == 6)])]
          if (length(indices_r) > 0) {
            symptomatic[indices_r] = round(triangle::rtriangle(
              n = length(indices_r),
              a = 5,
              b = 15,
              c = 10
            ), 0)
            ########## update recurrent indices
            update_2p = indices_r[which(indices_r %in% recur_1_idx)]
            update_1 = indices_r[which(!(
              indices_r %in% recur_2p_idx | indices_r %in% recur_1_idx
            ))]
            ## add patients with 2nd recur
            recur_2p_idx = c(recur_2p_idx, update_2p)
            ## remove patients with 1st recur who now have second
            recur_1_idx = recur_1_idx[!recur_1_idx %in% update_2p]
            ## add patients with first recur
            recur_1_idx = c(recur_1_idx, update_1)
            ## update results
            results$tot_symp_bc_recur = results$tot_symp_bc_recur + length(indices_r)
          }
          ## tran_7: from clinically resolved to susceptible
          indices_sus = clin_res_end[clin_res_end %in% (after_clin_res_location@i + 1)[which(after_clin_res_location@x == 7)]]
          if (length(indices_sus) > 0) {
            ## remove those going to sus from the recurrent tracker vectors
            recur_1_idx = recur_1_idx[!recur_1_idx %in% indices_sus]
            recur_2p_idx = recur_2p_idx[!recur_2p_idx %in% indices_sus]
            ## move to susceptible
            susceptible[indices_sus] = 1L
          }
        }
        ## tran_8: from asymptomatic infection to susceptible
        if (length(asymptomatic_end) > 0) {
          susceptible[asymptomatic_end] = 1L
        }
      }

      ## step 3: draw LOS for new symp ###############################################
      ## draw new LOS for patients who moved to symptomatic infection (new symp OR recurr.)
      # incubation_end is the indices of patients who are just moving to SI
      # indiced_r is the indices of patients who are moving to recurrent SI
      if (length(incubation_end) > 0 & (length(clin_res_end) > 0)) {
        new_los <- rep(NA, length(incubation_end) + length(indices_r))
        idx_now_in_symp = c(incubation_end, indices_r)
      } else if (length(clin_res_end) > 0) {
        ## and length(indices_r) > 0
        new_los <- rep(NA, length(indices_r))
        idx_now_in_symp = indices_r
      } else if (length(incubation_end) > 0) {
        new_los <- rep(NA, length(incubation_end))
        idx_now_in_symp = incubation_end
      } else {
        idx_now_in_symp <- numeric(0)
      }
      ## indices of new symp. patients not currently in a facility
      idx_symp_pat_nihosp = idx_now_in_symp[!idx_now_in_symp %in% (room_list$occup@x)]

      if (length(idx_now_in_symp) > 0) {
        ## new symptomatic who are in the hospital
        idx_symp_pat_in_hosp = idx_now_in_symp[idx_now_in_symp %in% (room_list$occup@x)]
        # length(idx_symp_pat_nihosp)
        # length(idx_now_in_symp)

        # draw los using CDI los function
        for (s in 1:length(idx_now_in_symp)) {
          # s = 1
          if (idx_now_in_symp[s] %in% idx_symp_pat_nihosp) {
            ## pat not in the hospital
            new_los[s] <-
              draw_CDI_los(
                current_day = 0,
                transfer = 1,
                md_cat = 6,
                ## mdc of prim dx CDI
                hcup_los = 0,
                hosp_id = pat_list$facility[idx_now_in_symp[s]]
              )
          } else if (idx_now_in_symp[s] %in% idx_symp_pat_in_hosp) {
            ## pat in the hospital
            # if(d == 57) {
            #   print(paste0("ADJ inhosp LOS. current_day: ", pat_list$los_sim[idx_now_in_symp[s]] - pat_list$days_rem[idx_now_in_symp[s]],
            #                " transfer: ", pat_list$tran_stat[idx_now_in_symp[s]],
            #                " md_cat: ", pat_list$mdc[idx_now_in_symp[s]],
            #                " hcup_los: ", pat_list$los_sim[idx_now_in_symp[s]],
            #                " hosp_id: ", pat_list$facility[idx_now_in_symp[s]]))
            # }
            new_los[s] <-
              draw_CDI_los(
                current_day = pat_list$los_sim[idx_now_in_symp[s]] - pat_list$days_rem[idx_now_in_symp[s]],
                transfer = pat_list$tran_stat[idx_now_in_symp[s]],
                md_cat = pat_list$mdc[idx_now_in_symp[s]],
                hcup_los = pat_list$los_sim[idx_now_in_symp[s]],
                hosp_id = pat_list$facility[idx_now_in_symp[s]]
              )

          }
        }
        ## update LOS details if already in the hospital
        ## indices of the vector `idx_now_in_symp` who are in the hospital
        idx_temp = which(idx_now_in_symp %in% idx_symp_pat_in_hosp)
        ### update LOS if in the hosp
        pat_list$los_sim[(idx_now_in_symp[idx_temp])] = new_los[idx_temp]
        ### update LOS rem if in the hosp
        pat_list$days_rem[(idx_now_in_symp[idx_temp])] = new_los[idx_temp]
      }
      if (length(idx_symp_pat_nihosp) > 0) {
        new_symp_pat = tibble(
          patid = idx_symp_pat_nihosp,
          hcup_id = pat_list$facility[idx_symp_pat_nihosp],
          adrgriskmortality = 1,
          los_sim = new_los[which(idx_now_in_symp %in% idx_symp_pat_nihosp)]
        )
      } else {
        new_symp_pat = tibble(
          patid = numeric(),
          hcup_id = numeric(),
          adrgriskmortality = numeric()
        )
      }

      ## step 4: pat days left in stay ###############################################
      ## update patient days left in stay
      # check patient discharges
      idx_to_discharge = idx_to_discharge[!idx_to_discharge %in% incubation_end]
      # `idx_to_discharge` obtained in step 10 (or initialization)

      # # check patient discharges
      #
      # ## indices where someone just completed last day
      # idx_to_discharge = (pat_list$days_rem@i + 1)[which(pat_list$days_rem@x == 1)]  # changed from == 1
      # ## update days left in stay for non discharged patients
      # idx_not_discharge = (pat_list$days_rem@i + 1)[pat_list$days_rem@x > 1]
      # pat_list$days_rem@x = pat_list$days_rem@x - 1

      ## indices where someone is being discharged & they are a transfer
      any_non_final_tran = (pat_list$tran_stat@i + 1)[which(pat_list$tran_stat@x == 3)]
      idx_to_tran = any_non_final_tran[any_non_final_tran %in% idx_to_discharge]
      ## discharge patients
      pat_list$days_rem[idx_to_discharge] = 0L

      ## step 5: discharge pat #######################################################
      ## discharged patients at the end of their LOS
      ## remove rest of dynamic vars from discharged pat
      pat_list$los_sim[idx_to_discharge] = 0L
      # pat_list$tran_stat[idx_to_discharge] = 0L
      pat_list$mdc[idx_to_discharge] = 0L
      pat_list$room_dwell[idx_to_discharge] = 0L
      # pat_list$facility[idx_to_discharge] = 0L    # for now, not setting this to zero in order to have facility if viz is constructed for symp CDI
      pat_list$risk_mor[idx_to_discharge] = 0L
      pat_list$new_symp_viz[idx_to_discharge] = 0L

      ## remove discharged patients from their rooms
      # room indices where someone is being discharged
      idx_discharge_rm = (room_list$occup@i + 1)[(room_list$occup@x %in% idx_to_discharge)]
      if (length(idx_to_discharge) != length(idx_discharge_rm)) {
        stop(
          "step 5 error: number of discharged patients is not equal to the number of discharged beds"
        )
      }
      # empty rooms
      room_list$occup[idx_discharge_rm] = 0L
      # finish patient discharge
      pat_list$viz_key[idx_to_discharge] = 0L

      ## step 6: room cleaning #######################################################
      ## rooms cleaning (terminal cleaning for discharges & daily cleaning
      ## for rest). Updates room contamination status if applicable
      # room indices where rooms are contam
      idx_contam = (room_list$contam@i + 1)[which(room_list$contam@x == 1)]
      if (length(idx_contam) > 0) {
        # room indices where someone is being discharged ( moved a few lines up)
        # idx_discharge_rm = (room_list$occup@i + 1)[which(room_list$occup@x %in% idx_to_discharge)]
        # rooms where patient was discharged and the room is contaminated
        idx_term_clean = idx_discharge_rm[idx_discharge_rm %in% idx_contam]
        # reset room contam with terminal clean prob
        room_list$contam[idx_term_clean] =
          as.integer(rbinom(
            n = length(idx_term_clean),
            size = 1,
            prob = (1 - 0.89)
          ))
        # daily cleaning for all other contaminated rooms where patient is still present
        idx_daily = idx_contam[!idx_contam %in% idx_discharge_rm]
        # reset room contam with daily clean prob
        room_list$contam[idx_daily] =
          as.integer(rbinom(
            n = length(idx_daily),
            size = 1,
            prob = (1 - 0.22)
          ))
      }

      ## step 7.1: admit patients ######################################################
      ## find all brand new patients with a visit
      ## new hcup patients with at least 1 visit already (inter-viz time is up)
      incoming_pat = (pat_list$days_to_next_viz@i + 1)[(pat_list$days_to_next_viz@x == 1)]
      incoming_pat = incoming_pat[!incoming_pat %in% room_list$occup@x]
      # same as: incoming_pat[incoming_pat %in% idx_to_discharge]

      # 7a: get new symptomatic patients not in hospital to admit
      # (df created at step 3: new_symp_pat)
      symp_pat_to_enter_early_ortoday = new_symp_pat$patid[pat_list$tran_stat[new_symp_pat$patid] == 3]
      ## remove patients from coming in early if already slated to come in today
      # symp_pat_to_enter_early = symp_pat_to_enter_early[!symp_pat_to_enter_early %in% incoming_pat]

      df_enter_early = tibble(
        patid = symp_pat_to_enter_early_ortoday,
        viz_num = pat_list$viz_num[patid] + 1
      )
      ## take early enter patients out of new symp, they enter with transfers (even tho they are newly symp)
      new_symp_pat = new_symp_pat |>
        filter(!patid %in% symp_pat_to_enter_early_ortoday) |>
        filter(!patid %in% queue_viz_keys_df$patid)
      ## adjust vector of symp patients
      remove_idx = c(symp_pat_to_enter_early_ortoday, queue_viz_keys_df$patid)
      idx_symp_pat_nihosp = idx_symp_pat_nihosp[!idx_symp_pat_nihosp %in% remove_idx]
      ## get secondary viz patients (including new symp coming early)
      scndry_viz_pats = tibble(
        patid = incoming_pat,
        viz_num = pat_list$viz_num[patid] + 1
      ) |>
        filter(!patid %in% symp_pat_to_enter_early_ortoday) |>
        bind_rows(df_enter_early) |>
        inner_join(hcup, join_by(patid, viz_num))
      # pat_w_next = pat_list$days_to_next_viz@i + 1
      # pat_w_next_nihosp = pat_w_next[!pat_w_next %in% (room_list$occup@x)]
      # pat_w_next_nihosp[pat_w_next_nihosp %in% which(pat_list$days_to_next_viz == 1)]

      ## admit patients to HCFs & assign rooms
      # 7b: get transfer patients to admit (from hcup & currently in the hospital)
      # get keys of the discharged patients with a transfer next
      viz_key_next_tran_seg <-
        scndry_viz_pats |>
        filter(subseq_tran_seg == 1) |>
        select(patid, viz_key, hcup_id, adrgriskmortality)

      # 7c: get queue patients to admit
      # queue_viz_keys_df
      ## adjust LOS for symp queue patients
      # browser()
      adj_queue_los_df = queue_viz_keys_df |>
        select(-los_sim) |>
        filter(patid %in% new_symp_pat$patid) |>
        inner_join(hcup, join_by(patid, viz_key, hcup_id, adrgriskmortality)) |>
        select(patid, viz_key, hcup_id, los_sim, tran_seg, mdc)
      if (nrow(adj_queue_los_df) > 0) {
        new_los <- rep(NA, nrow(adj_queue_los_df))
        for (s in 1:nrow(adj_queue_los_df)) {
          # if(d == 57) {
          #   print(paste0("ADJ Q LOS. current_day: ", pat_list$los_sim[idx_now_in_symp[s]] - pat_list$days_rem[idx_now_in_symp[s]],
          #                " transfer: ", pat_list$tran_stat[idx_now_in_symp[s]],
          #                " md_cat: ", pat_list$mdc[idx_now_in_symp[s]],
          #                " hcup_los: ", pat_list$los_sim[idx_now_in_symp[s]],
          #                " hosp_id: ", pat_list$facility[idx_now_in_symp[s]]))
          # }
          new_los[s] = draw_CDI_los(
            current_day = 0,
            transfer = adj_queue_los_df$tran_seg[s],
            md_cat = adj_queue_los_df$mdc[s],
            hcup_los = adj_queue_los_df$los_sim[s],
            hosp_id = adj_queue_los_df$hcup_id[s]
          )
        }
        adj_queue_los_df = adj_queue_los_df |> mutate(adj_los = new_los)
      }
      # 7d.1: get available rooms per facility
      update_available_rooms(rm_list = room_list, hcupids_vec = hcup_id_vec)
      ## R func that updates environ automatically, list of vectors of avail rooms for each facility

      ## check if any patients are dead
      idx_dead = (death@i + 1)

      to_admit_pat_7_1 =
        viz_key_next_tran_seg |>       # from 7b
        bind_rows(queue_viz_keys_df) |>   # 7c: join_by(key, ahaid, adrgriskmortality)
        bind_rows(new_symp_pat) |>        # 7b
        filter(!patid %in% idx_dead)

      ## check if any patients to be admitted are actually already in the hospital
      idx_already_inhosp = to_admit_pat_7_1$patid[to_admit_pat_7_1$patid %in% room_list$occup@x]

      if (length(idx_already_inhosp) > 0) {
        print(
          "trying to admit someone that is currently in the hosp->keeping them in queue, l566"
        )
        send_to_q_df = to_admit_pat_7_1 |>
          filter(patid %in% idx_already_inhosp) |>
          select(patid, viz_key, hcup_id, adrgriskmortality)
        ## put in queue
        queue_viz_keys_df = bind_rows(queue_viz_keys_df, send_to_q_df)
        ## adjust patients to admit
        to_admit_pat_7_1 =
          to_admit_pat_7_1 |>
          filter(!viz_key %in% send_to_q_df$viz_key)
      }

      ## num of patients to admit
      to_admit_summ =
        to_admit_pat_7_1 |>
        count(hcup_id, adrgriskmortality) |> ## num pat per facility in each risk group
        bind_rows(tibble(adrgriskmortality = 1:4, n = rep(0, 4))) |>
        pivot_wider(
          names_from = adrgriskmortality,
          names_prefix = "risk_",
          values_from = n
        ) |>
        ungroup() |>
        replace_na(list(risk_1 = 0, risk_2 = 0, risk_3 = 0, risk_4 = 0)) |>
        filter(!is.na(hcup_id)) |>
        select(hcup_id, risk_1, risk_2, risk_3, risk_4) |>
        left_join(
          tibble(
            fid = hcup_id_vec,
            n_icu = sapply(icu_avail, length),
            n_non = sapply(non_avail, length)
          ),
          join_by(hcup_id == fid)
        ) |>
        data.frame()

      n_in_queue =
        num_in_q_by_risk(risks_beds = to_admit_summ) ## list with num by risk / facility
      ## will everyone fit?
      ppl_to_queue = sapply(n_in_queue, any) |> any()
      leave_in_queue = c() ## initalize vector to hold queue patient visits
      # print(paste0("total patient to queue: ", sum(sapply(n_in_queue, sum))))
      if (ppl_to_queue == TRUE) {
        print("patients going to queue at queue/tran/symp patients, l603")
        for (f in 1:length(n_in_queue)) {
          # n ppl to be sent to queue for facility f
          n_to_q_for_f = sum(n_in_queue[[f]])
          if (n_to_q_for_f > 0) {
            ## people  to be sent to queue
            n_in_prev_q = sum(queue_viz_keys_df$hcup_id == names(n_in_queue)[f])
            n_symp_in_q = sum(adj_queue_los_df$hcup_id == names(n_in_queue)[f])
            n_in_prev_q = n_in_prev_q -  n_symp_in_q
            if (n_in_prev_q >= n_to_q_for_f) {
              ## if enough to just leave queue people in queue
              # n_to_leave = n_in_prev_q - n_to_q_for_f
              key_to_leave_in_q =
                queue_viz_keys_df |>
                filter(hcup_id == names(n_in_queue)[f]) |>
                filter(!patid %in% adj_queue_los_df$patid) |> ## so symp queue patients enter
                filter(patid %in% to_admit_pat_7_1$patid) |>
                slice_sample(n = n_to_q_for_f) |>
                pull(patid)
              leave_in_queue = c(leave_in_queue, key_to_leave_in_q)
            } else {
              ## need new queue patients from old queue and new symp patients
              q_and_symp_pat =
                to_admit_pat_7_1 |>
                filter(hcup_id == names(n_in_queue)[f]) |>
                filter(!patid %in% adj_queue_los_df$patid) |>
                mutate(
                  viz_cat = case_when(
                    !is.na(los_sim) ~ "symp",
                    patid %in% queue_viz_keys_df$patid ~ "queue",
                    TRUE ~ "tran"
                  )) |>
                filter(viz_cat != "tran")
              # print(q_and_symp_pat)
              tot_q_and_symp_pat = nrow(q_and_symp_pat)
              if(tot_q_and_symp_pat >= n_to_q_for_f) {
                ## select queue
                df_queue = q_and_symp_pat |>
                  slice_sample(n = n_to_q_for_f)
                key_to_leave_in_q = pull(df_queue, patid)
                # print(key_to_leave_in_q)
                leave_in_queue = c(leave_in_queue, key_to_leave_in_q)
                ## adjust `idx_symp_pat_nihosp`
                symp_adj = df_queue |> filter(viz_cat == "symp") |> pull(patid)
                idx_symp_pat_nihosp = idx_symp_pat_nihosp[!idx_symp_pat_nihosp %in% symp_adj]
                ## not adjusting `symptomatic` days status and letting them enter a day or so later
              } else {
                stop("error at step 7c: more people from transfers need to go in the queue!!!!!")
            }
          }
        }
        }
      }
      # if(d == 42) {browser()}
      # print(paste0("dim queue before: ", dim(queue_viz_keys_df)))
      ## adjust queue. `los_sim` var should only be included in queue patients
          ## if they are a new symp viz, otherwise will mess with when/how those patients are admitted
      queue_viz_keys_df =
        queue_viz_keys_df |>
        filter(patid %in% leave_in_queue)
      # print(paste0("dim queue after: ", dim(queue_viz_keys_df)))
      ## adjust patients to admit
      # print(paste0("dim admit pat df before: ", dim(to_admit_pat_7_1)))
      to_admit_pat_7_1 =
        to_admit_pat_7_1 |>
        filter(!patid %in% leave_in_queue)
      # print(paste0("dim admit pat df after: ", dim(to_admit_pat_7_1)))
      # 7d.3: assign rooms for patients from 7a, 7b, 7c
      new_rooms <- vector("list", length(hcup_id_vec))
      # assign rooms
      for (i in 1:length(hcup_id_vec)) {
        ## this for loop is fast
        temp_to_admit =
          to_admit_pat_7_1 |>
          filter(hcup_id == hcup_id_vec[i]) |>
          arrange(desc(adrgriskmortality)) |>
          select(patid, adrgriskmortality)
        current_icu = icu_avail[[as.character(hcup_id_vec[i])]]
        current_non = non_avail[[as.character(hcup_id_vec[i])]]
        if(nrow(temp_to_admit) > length(current_icu) + length(current_non)){
          print(paste0("i=", i, "; tot pat: ", nrow(temp_to_admit), "; tot rooms:", length(current_icu) + length(current_non),"; icu: ", length(current_icu), "; non: ", length(current_non)))
        }
        new_rooms[[i]] =
          rcpp_assign_rooms_cpp_seed(
            pat_risks = temp_to_admit,
            icu = current_icu,
            non = current_non,
            seed = SEED
          )
      }
      new_rooms <- bind_rows(new_rooms)
      # print(dim(new_rooms))
      ## check all the right people: all(sort(to_admit_pat_7_1$patid) == sort(new_rooms$patid))
      # update occupancy of rooms from 7d.3 assignments
      room_list$occup[new_rooms$assigned_room] = new_rooms$patid
      ## adjust LOS if queue pat is now SYMP
      if(nrow(adj_queue_los_df) > 0) {
        to_admit_df <- hcup |>
          filter(viz_key %in% to_admit_pat_7_1$viz_key & !viz_key %in% all_admit_keys) |>
          left_join(adj_queue_los_df |> select(patid, viz_key, los_sim, adj_los)) |>
          mutate(los_sim = if_else(!is.na(adj_los), adj_los, los_sim)) |>
          select(-adj_los)
        print(paste0("There are", nrow(adj_queue_los_df), " symp pat from the queue on day ", d))
      } else {
        to_admit_df <- hcup |>
          filter(viz_key %in% to_admit_pat_7_1$viz_key & !viz_key %in% all_admit_keys)
      }
      # add patient information for transfer/queue admits
      all_admit_keys = c(all_admit_keys, to_admit_df$viz_key) ## dynamic vector of all used keys
      to_admit_pat = to_admit_df$patid
      pat_list$los_sim[to_admit_pat] = to_admit_df$los_sim
      pat_list$days_rem[to_admit_pat] = to_admit_df$los_sim
      pat_list$tran_stat[to_admit_pat] = to_admit_df$tran_seg
      pat_list$tran_sub_fac_diff[to_admit_pat] = to_admit_df$fac_diff_flag
      pat_list$tran_days_bn[to_admit_pat] = to_admit_df$days2next
      pat_list$mdc[to_admit_pat] = to_admit_df$mdc
      pat_list$facility[to_admit_pat] = to_admit_df$hcup_id
      pat_list$risk_mor[to_admit_pat] = to_admit_df$adrgriskmortality
      pat_list$viz_key[to_admit_pat] = to_admit_df$viz_key
      pat_list$viz_num[to_admit_pat] = to_admit_df$viz_num
      ## add days to next viz, if not NA
      next_viz_dat = to_admit_df |>
        filter(!is.na(days2next))
      pat_list$days_to_next_viz[next_viz_dat$patid] = next_viz_dat$days2next
      ## add days since prev viz for revisits, if not NA
      next_viz_dat = to_admit_df |>
              filter(!is.na(revisit_days_since))
      pat_list$revisit_days_since[next_viz_dat$patid] = next_viz_dat$revisit_days_since

      # add patient info for newly symp pat admits  (Condition used to be
      # length(idx_symp_pat_nihosp) > 0 will need. Was adjusted here b/c queue
      # needed to include new symp not in hosp)
      if (sum(!is.na(to_admit_pat_7_1$los_sim)) > 0) {
        new_symp_admit = to_admit_pat_7_1 |> filter(!is.na(los_sim))
        if(any(!is.na(new_symp_admit$viz_key))) {stop("error step 7d: viz key patients are trying to be admitted with the new symptomatic visits")}
        idx_new_symp_admit = new_symp_admit$patid
        ## this is updated in step 3 for incubation period
        pat_list$los_sim[idx_new_symp_admit] = new_symp_admit$los_sim
        pat_list$days_rem[idx_new_symp_admit] = new_symp_admit$los_sim
        pat_list$tran_stat[idx_new_symp_admit] = 1L
        pat_list$mdc[idx_new_symp_admit] = 6L
        # pat_list$facility     # this stays as is from previous visit
        pat_list$risk_mor[idx_new_symp_admit] = 1L # could assign risk differently but
        # ~90% pat are in risk 1-3 so doesn't make a big difference & only
        # using this for assigning icu/non at admission
        pat_list$new_symp_viz[idx_new_symp_admit] = 1L
      }

      ## update disease status for any first time visit patients (from queue)
      if (any(to_admit_df$viz_num == 1)) {
        first_visits = filter(to_admit_df, viz_num == 1)
        ## patients with CDI status
        cdi_admits = first_visits |>
          select(patid) |>
          left_join(uniq_pat_dat |> select(patid, class, class_sim),
                    join_by(patid)) |>
          filter(class_sim == "CDI")
        ## patients entering at symptomatic
        idx_to_symp = cdi_admits |>
          filter(class %in% c("R", "HAI", "T")) |>
          pull(patid)
        if (length(idx_to_symp) > 0) {
          symptomatic[idx_to_symp] =
            round(triangle::rtriangle(n = length(idx_to_symp), a = 5, b = 15, c = 10), 0)
        }
        ## patients entering at incubation
        idx_to_incub = cdi_admits |>
          filter(class == "CAI") |>
          pull(patid)
        if (length(idx_to_incub) > 0) {
          incubation[idx_to_incub] =
            round(triangle::rtriangle(n = length(idx_to_incub), a = 1, b = 10, c = 4))
        }
        ## non CDI patients
        ## set probability of noCDI patients to be asymp at first viz
        maybe_asymp_pat =
          first_visits |>
          filter(viz_num == 1 & !patid %in% cdi_admits$patid) |>
          pull(patid)
        idx_pat_asymp =
          maybe_asymp_pat[which(rbinom(length(maybe_asymp_pat), size = 1, prob = 0.00785) == 1)]
        asymptomatic[idx_pat_asymp] = round(rweibull(length(idx_pat_asymp), 1, 43), 0) + 1
        idx_pat_sus = maybe_asymp_pat[!maybe_asymp_pat %in% idx_pat_asymp]
        susceptible[idx_pat_sus] = 1L
      }

      ## 7e: ADMIT PART 2: admit new patients from HCUP #######################
      # 7e.1: new hcup pat - get available rooms
      update_available_rooms(rm_list = room_list, hcupids_vec = hcup_id_vec)  ## updates envir automatically

      ## first time patients
      first_viz_pats = uniq_pat_dat |>
        filter(aday_sim == d) |>
        select(patid) |>
        mutate(viz_num = 1) |>
        inner_join(hcup, join_by(patid, viz_num))

      ## secondary non-tran patients
      sec_viz_nont_pats = scndry_viz_pats |>
        filter(subseq_tran_seg == 0 &
                 !patid %in% idx_dead & !viz_key %in% all_admit_keys) #|>
      # select(patid, viz_key, hcup_id, adrgriskmortality)

      # 7e.2: new hcup pat - will everyone fit? to queue?
      to_admit_df <- bind_rows(first_viz_pats, sec_viz_nont_pats)
      ## check if any patients to be admitted are actually already in the hospital
      idx_already_inhosp = to_admit_df$patid[which(to_admit_df$patid %in% room_list$occup@x)]
      if (length(idx_already_inhosp) > 0) {
        print(
          "trying to admit someone that is currently in the hosp->sending them to queue, l608"
        )
        send_to_q_df = to_admit_df |>
          filter(patid %in% idx_already_inhosp) |>
          select(patid, viz_key, hcup_id, adrgriskmortality)
        ## put in queue
        queue_viz_keys_df = bind_rows(queue_viz_keys_df, send_to_q_df)
        ## adjust patients to admit
        to_admit_df =
          to_admit_df |>
          filter(!viz_key %in% send_to_q_df$viz_key)
      }


      ## num of patients to admit & n rooms avail by icu/non
      to_admit_summ =
        to_admit_df |>
        count(hcup_id, adrgriskmortality) |> ## num pat per facility in each risk group
        bind_rows(tibble(adrgriskmortality = 1:4, n = rep(0, 4))) |>
        pivot_wider(
          names_from = adrgriskmortality,
          names_prefix = "risk_",
          values_from = n
        ) |>
        ungroup() |>
        replace_na(list(risk_1 = 0, risk_2 = 0, risk_3 = 0, risk_4 = 0)) |>
        filter(!is.na(hcup_id)) |>
        select(hcup_id, risk_1, risk_2, risk_3, risk_4) |>
        left_join(
          tibble(
            fid = hcup_id_vec,
            n_icu = sapply(icu_avail, length),
            n_non = sapply(non_avail, length)
          ),
          join_by(hcup_id == fid)
        ) |>
        data.frame()

      ## will everyone fit?
      n_in_queue =
        num_in_q_by_risk(risks_beds = to_admit_summ) ## list with num by risk / facility
      ppl_to_queue = sapply(n_in_queue, any) |> any()
      leave_in_queue = c()
      if (ppl_to_queue == TRUE) {
        # fac_w_q = names(which(fac_w_q == TRUE))
        print("patients going to queue at new hcup patients, l637")
        for (f in 1:length(n_in_queue)) {
          # n ppl to be sent to queue for facility f
          n_to_q_for_f = sum(n_in_queue[[f]])
          # send_to_q_df = tibble(patid = numeric(), viz_key = numeric(), hcup_id = numeric(), adrgriskmortality = numeric())
          if (n_to_q_for_f > 0) {
            ## people  to be sent to queue
            send_to_q_df =
              to_admit_df |>
              filter(hcup_id == names(n_in_queue)[f]) |>
              slice_sample(n = n_to_q_for_f) |>
              select(patid, viz_key, hcup_id, adrgriskmortality)
            # adjust the queue
            queue_viz_keys_df = bind_rows(queue_viz_keys_df, send_to_q_df)
          }
        }
      }
      ## adjust patients to admit for queue (remove patients sent to queue)
      to_admit_df =
        to_admit_df |>
        filter(!viz_key %in% queue_viz_keys_df$viz_key)
      # 7e.3: new hcup pat - assign rooms & admit patients
      # assign rooms
      new_rooms <- vector("list", length(hcup_id_vec))
      for (i in 1:length(hcup_id_vec)) {
        ## this for loop is fast
        temp_to_admit =
          to_admit_df |>
          filter(hcup_id == hcup_id_vec[i]) |>
          arrange(desc(adrgriskmortality)) |>
          select(patid, adrgriskmortality) |>
          mutate(adrgriskmortality = if_else(adrgriskmortality == 0 | is.na(adrgriskmortality), 1, adrgriskmortality))
        current_icu = icu_avail[[as.character(hcup_id_vec[i])]]
        current_non = non_avail[[as.character(hcup_id_vec[i])]]
        if(nrow(temp_to_admit) > length(current_icu) + length(current_non)){
          print(paste0("i=", i, "; tot pat: ", nrow(temp_to_admit), "; tot rooms:", length(current_icu) + length(current_non),"; icu: ", length(current_icu), "; non: ", length(current_non)))
        }
        new_rooms[[i]] =
          rcpp_assign_rooms_cpp_seed(
            pat_risks = temp_to_admit,
            icu = current_icu,
            non = current_non,
            seed = SEED
          )
      }
      new_rooms <- bind_rows(new_rooms)
      # set rooms as occupied
      room_list$occup[new_rooms$assigned_room] = new_rooms$patid
      # # add patient information for new hcup admits
      all_admit_keys = c(all_admit_keys, to_admit_df$viz_key) ## dynamic vector of all used keys
      # admit patients
      to_admit_pat = to_admit_df$patid
      pat_list$los_sim[to_admit_pat] = to_admit_df$los_sim
      pat_list$days_rem[to_admit_pat] = to_admit_df$los_sim
      pat_list$tran_stat[to_admit_pat] = to_admit_df$tran_seg
      pat_list$tran_sub_fac_diff[to_admit_pat] = to_admit_df$fac_diff_flag
      pat_list$tran_days_bn[to_admit_pat] = to_admit_df$days2next
      pat_list$mdc[to_admit_pat] = to_admit_df$mdc
      pat_list$facility[to_admit_pat] = to_admit_df$hcup_id
      pat_list$risk_mor[to_admit_pat] = to_admit_df$adrgriskmortality
      pat_list$viz_key[to_admit_pat] = to_admit_df$viz_key
      pat_list$viz_num[to_admit_pat] = to_admit_df$viz_num
      ## add days to next viz, if not NA
      next_viz_dat = to_admit_df |>
        filter(!is.na(days2next))
      pat_list$days_to_next_viz[next_viz_dat$patid] = next_viz_dat$days2next
      ## add days since prev viz for revisits, if not NA
      next_viz_dat = to_admit_df |>
        filter(!is.na(revisit_days_since))
      pat_list$revisit_days_since[next_viz_dat$patid] = next_viz_dat$revisit_days_since

      ## update disease status for first-visit patients (viz_num = 1)
      ## patients with viz_num > 1 do not get their disease status updated here
      ## patients with CDI status
      cdi_admits = to_admit_df |>
        filter(viz_num == 1) |>
        select(patid) |>
        left_join(uniq_pat_dat |> select(patid, class, class_sim),
                  join_by(patid)) |>
        filter(class_sim == "CDI")
      ## patients entering at symptomatic
      idx_to_symp = cdi_admits |>
        filter(class %in% c("R", "HAI", "T")) |>
        pull(patid)
      if (length(idx_to_symp) > 0) {
        symptomatic[idx_to_symp] = round(triangle::rtriangle(
          n = length(idx_to_symp),
          a = 5,
          b = 15,
          c = 10
        ), 0)
      }
      ## patients entering at incubation
      idx_to_incub = cdi_admits |>
        filter(class == "CAI") |>
        pull(patid)
      if (length(idx_to_incub) > 0) {
        incubation[idx_to_incub] = round(triangle::rtriangle(
          n = length(idx_to_incub),
          a = 1,
          b = 10,
          c = 4
        ))
      }
      ## non CDI patients
      ## set probability of noCDI patients to be asymp at first viz
      maybe_asymp_pat =
        to_admit_df |>
        filter(viz_num == 1 & !patid %in% cdi_admits$patid) |>
        pull(patid)
      idx_pat_asymp = maybe_asymp_pat[which(rbinom(
        length(maybe_asymp_pat),
        size = 1,
        prob = 0.00785
      ) == 1)]
      asymptomatic[idx_pat_asymp] = round(rweibull(length(idx_pat_asymp), 1, 43), 0) + 1
      idx_pat_sus = maybe_asymp_pat[!maybe_asymp_pat %in% idx_pat_asymp]
      susceptible[idx_pat_sus] = 1L

      ## find symptomatic patients who are recurrent
      # recur_pat = after_clin_res_location@i[which(after_clin_res_location@x == 6)] + 1
      recur_pat = c(recur_1_idx, recur_2p_idx)
      recur_pat_today = recur_pat[recur_pat %in% (symptomatic@i + 1)]

      ## step 8: pat mvt patterns for day #############################################
      ## draw patient movement patterns for the day for all patients in hospitals

      # 8a: sample the number of mvts per patient for the day

      # indices of patients in the hospital
      pat_in_hosp_idx = (pat_list$days_rem@i + 1)    # don't use facility Matrix b/c not reset upon discharge
      # order of patients rooms in the hospit
      order_pat_rooms_idx = order(room_list$occup@x)
      # indices of occupied rooms in the hospital
      occup_rm_idx = room_list$occup@i + 1
      # icu flags for occupied rooms
      pat_room_type = as.integer((room_list$icu[occup_rm_idx])[order_pat_rooms_idx]) ## order by asc patient idx
      # sample the number of icu/non movements for each patient
      pat_mvt_num_draws = sample_day_mvts_cpp_seed(
        los = pat_list$los_sim[pat_in_hosp_idx],
        cur_room_type = pat_room_type,
        seed = SEED
      )
      # if(d == 44) {browser()}
      if (length(pat_in_hosp_idx) != length(order_pat_rooms_idx)) {

        tibble(occup = room_list$occup@x) |> count(occup, sort = TRUE) |> print(n = 10)
        tibble(pat = pat_in_hosp_idx) |> count(pat, sort = TRUE) |> print(n = 10)
        # who doesn't have a room??
        pat_in_hosp_idx[!pat_in_hosp_idx %in% room_list$occup@x] |> print()
        ## patients in a room but not in the "patients in hosp index" (should have been discharged)
        room_list$occup@x[!room_list$occup@x %in% pat_in_hosp_idx] |> print()
        stop(
          "step 8a error: number of hospitalizated patients is not equal to the number of occupied beds"
        )
      }
      ## investigate:
      # tibble(occup = room_list$occup@x) |> count(occup, sort = TRUE)
      # tibble(pat = pat_in_hosp_idx) |> count(pat, sort = TRUE)
      # # who doesn't have a room??
      pat_in_hosp_idx[!pat_in_hosp_idx %in% room_list$occup@x]
      # ## patients in a room but not in the "patients in hosp index" (should have been discharged)
      # room_list$occup@x[!room_list$occup@x %in% pat_in_hosp_idx]

      # 8b: create data frame of patient movements
      # columns: id, room_type(0/1), time_block(1:4); one row per patient mvt in the day
      pat_mvts_day_df <- create_pat_mov_df(ids = pat_in_hosp_idx,
                                           start_loc = pat_room_type,
                                           drawn_moves = pat_mvt_num_draws)
      # output col: patid, room_type, time_block
      ## step 9: run four time blocks ################################################
      ## run 4 time blocks for the day (6 hours each)
      for (tb in 1:4) {
        # tb = 4
        ## step 9a: patients move rooms if applicable #######################################
        ## patients to move
        pat_to_move = pat_mvts_day_df |>
          filter(time_block == tb) |>
          select(patid, room_type)
        ## add hcup facility to get room from
        pat_to_move = pat_to_move |>
          mutate(hcup_id = pat_list$facility[pat_to_move$patid])
        ## reset room occupancy of moving patients
        # room indices where someone is being moved
        idx_move_rm = (room_list$occup@i + 1)[which(room_list$occup@x %in% pat_to_move$patid)]
        room_list$occup[idx_move_rm] = 0L
        ## update available rooms
        update_available_rooms(rm_list = room_list, hcupids_vec = hcup_id_vec)
        ## draw new rooms for moving patients
        fac_w_move = unique(pat_to_move$hcup_id)
        new_rooms <- vector("list", length(fac_w_move))
        # new_rooms <- vector("list", length(hcup_id_vec))
        for (f in 1:length(fac_w_move)) {
          ## this for loop is fast
          # f = 1
          new_rooms[[f]] <-
            pat_to_move |>
            filter(hcup_id == fac_w_move[f]) |>
            arrange(desc(room_type)) |>
            select(patid, room_type) |>
            move_rooms_cpp_seed(
              pat_rm_type = _,
              icu = icu_avail[[as.character(fac_w_move[f])]],
              non = non_avail[[as.character(fac_w_move[f])]],
              seed = SEED
            )
        }
        new_rooms <- bind_rows(new_rooms)
        ## set new rooms as occupied
        room_list$occup[new_rooms$assigned_room] = new_rooms$patid

        ## step 9a.1: terminal cleaning of rooms where patients moved out of
        ## rooms that may need to be cleaned: idx_move_rm
        # room indices for any room that is contam
        idx_contam = (room_list$contam@i + 1)[which(room_list$contam@x == 1)]
        if (length(idx_contam) > 0) {
          # rooms where patient was moved and the room is contaminated
          idx_term_clean = idx_move_rm[idx_move_rm %in% idx_contam]
          # reset room contam with terminal clean prob
          room_list$contam[idx_term_clean] =
            as.integer(rbinom(
              n = length(idx_term_clean),
              size = 1,
              prob = (1 - 0.89)
            ))
        }

        ## step 9a.2: update room dwell time
        # pat_move_idx = pat_to_move$patid
        # pat_idx = (pat_list$days_rem@i + 1)
        # pat_stay_idx = pat_idx[!pat_idx %in% pat_move_idx]
        # pat_list$room_dwell[pat_move_idx] = 360L
        # pat_list$room_dwell[pat_stay_idx] = pat_list$room_dwell[pat_stay_idx] + 360L

        ## step 9b: create HCW mvt matrix for time block ###########################
        # d = 5; tb = 2
        hcw_mvts_block <-
          hcw_mvts |>
          filter(day_counter == d & time_block == tb) |>
          select(rid, hid, total_sec)
        ## dim: 5036    3
        tb_hcw_room_inter <- vector("list", length(hcup_id_vec))
        names(tb_hcw_room_inter) = hcup_id_vec
        for (f in 1:length(hcup_id_vec)) {
          # f = 1
          f_rooms = rooms |>
            filter(hcup_id == hcup_id_vec[f]) |>
            select(rid, rid_uniq)
          ## dim 88 x 4
          df <- hcw_mvts_block |>
            filter(rid %in% f_rooms$rid) |>
            rename(hid_ss = hid) |>
            left_join(f_rooms, join_by(rid)) |>
            left_join(hcws |> filter(hcup_id == hcup_id_vec[f]),
                      join_by(hid_ss == hid)) |>
            # left_join(hcw_lookup |> select(hid, prop_soap_in, prop_soap_out), join_by(hid_ss == hid)) |>
            select(
              # rid,
              rid_uniq,
              hid_uniq,
              total_sec,
              prop_soap_in,
              prop_soap_out
            )
          # df
          tb_hcw_room_inter[[f]] = df
        }
        tb_hcw_room_inter = bind_rows(tb_hcw_room_inter)
        ## steps 9c/d: update room contam probs ########################################

        ## step 9c: update room comtam prob on disease status of HCWs & patients

        # set up patient indicators based on research questions

        # RESEARCH Q1
        if (incl_col_pat_trans == 1) {
          # as-is/normal scenario
          idx_pat_col = c(incubation@i + 1, clin_res@i + 1, asymptomatic@i + 1)
        } else {
          # asymptomatic patients are dropped from colonized patient spread
          idx_pat_col = numeric()
        }

        if (incl_asymp_pat_trans == 1) {
          # as-is/normal scenario
          idx_pat_col = c(incubation@i + 1, clin_res@i + 1, asymptomatic@i + 1)
        } else {
          # asymptomatic patients are dropped from colonized patient spread
          idx_pat_col = c(incubation@i + 1, clin_res@i + 1)
        }

        if (incl_asymp_hcw_trans == 1) {
          # as-is/normal scenario
          idx_hcw_col = (hcw_list$colonized@i + 1)
        } else {
          # no hcws are colonized
          idx_hcw_col <- numeric()
        }

        # RESEARCH Q2
        if (incl_recur_trans == 1) {
          # as-is / normal
          idx_pat_symp = (symptomatic@i + 1)
        } else {
          # for symp transmission, only considering non-recurrent cases
          idx_pat_symp = (symptomatic@i + 1)[!(symptomatic@i + 1) %in% recur_pat_today]
        }

        # RESEARCH Q3
        # as-is/normal scenario: incl_transfer_pat_trans == 1 & incl_non_transfer_pat_transmission == 1
        # no change to idx_pat_symp & idx_pat_col vectors
        if (incl_transfer_pat_trans == 0) {
          # exclude transfer patient transmission, remove from idx vectors
          transfers = (pat_list$tran_stat@i[pat_list$tran_stat@x %in% c(2, 3)]) + 1
          # pat_list$tran_sub_fac_diff
          # pat_list$tran_days_bn
          idx_pat_symp = idx_pat_symp[!idx_pat_symp %in% transfers]
          idx_pat_col = idx_pat_col[!idx_pat_col %in% transfers]
        }

        # RESEARCH Q4
        if(incl_revisit_pat_trans == 0) {
          # rq_tran_days_bn
          revisits = which(pat_list$revisit_days_since <= rq_tran_days_bn)
          ## stop transmission from revisit patients
          idx_pat_symp = idx_pat_symp[!idx_pat_symp %in% revisits]
          idx_pat_col = idx_pat_col[!idx_pat_col %in% revisits]
        }

        # RESEARCH Q5
        if (incl_symp_trans == 1) {
          # as-is / normal
          idx_pat_symp = (symptomatic@i + 1)
        } else {
          # eliminate symp transmission
          idx_pat_symp = numeric()
        }

        #9c.1: prob_room_symp_pat: symptomatic infection
        # prob_room_symp_pat = 1 - exp(-pat_list$room_dwell[idx_pat_symp] * tau_symp_room)
        prob_room_symp_pat = rep(1 - exp(-360 * tau_symp_room), time = length(idx_pat_symp))
        #9c.2: prob_room_col_pat: incubation per, clin res, asymp inf
        # prob_room_col_pat = 1 - exp(-pat_list$room_dwell[idx_pat_col] * 0.59 * tau_symp_room)
        prob_room_col_pat = rep(1 - exp(-360 * 0.59 * tau_symp_room), times = length(idx_pat_col))
        ## room numbers of symp or colonized patients (need patients in correct room)
        rm_idx_check = c(idx_pat_symp, idx_pat_col) %in% room_list$occup@x
        # if(sum(rm_idx_check) != length(rm_idx_check)) {print("check if it is correct that not everyone symp/col pat has a room!!")}

        pat_room_prob_df =
          tibble(rid_uniq = (room_list$occup@i + 1)[which(room_list$occup@x %in% c(idx_pat_symp, idx_pat_col))],
                 patid = room_list$occup[rid_uniq]) |>
          left_join(tibble(
            patid = c(idx_pat_symp, idx_pat_col),
            prob_contam_fr_patdis = c(prob_room_symp_pat, prob_room_col_pat)
          ),
          join_by(patid)) |>
          select(-patid)

        #9c.3: prob_room_asymp_hcw (MORE COMPLEX b/c mult hcw per room )
        # idx_hcw_col = (hcw_list$colonized@i + 1) # move earlier with indicators
        idx_hcw_contam = (hcw_list$contam@i + 1)
        room_contam_probs =
          tb_hcw_room_inter |>
          mutate(
            prob_contam_col = if_else(
              hid_uniq %in% idx_hcw_col,
              1 - exp(-total_sec * 0.59 * tau_symp_room),
              NA
            ),
            prob_contam_contam = if_else(
              hid_uniq %in% idx_hcw_contam,
              (1 - prop_soap_in) * (1 - exp(-total_sec * tau_hcw_room)),
              NA
            )
          ) |>
          group_by(rid_uniq) |>
          summarise(
            prob_room_nocontam_hcw_col = prod(1 - prob_contam_col, na.rm = TRUE),
            prob_room_nocontam_hcw_contam = prod(1 - prob_contam_contam, na.rm = TRUE),
            .groups = "drop"
          ) |>
          left_join(pat_room_prob_df, join_by(rid_uniq)) |>
          replace_na(list(prob_contam_fr_patdis = 0)) |>
          mutate(
            prob_contam = 1 - ((1 - prob_contam_fr_patdis) * prob_room_nocontam_hcw_col * prob_room_nocontam_hcw_contam
            )
          ) |>
          # mutate(
          #   contam_fr_hcwdis = rbinom(n = length(prob_room_nocontam_hcw_col), size = 1, prob = (1 - prob_room_nocontam_hcw_col)),
          #   contam_fr_hcwcon = rbinom(n = length(prob_room_nocontam_hcw_contam), size = 1, prob = (1 - prob_room_nocontam_hcw_contam)),
          #   contam_fr_patdis = rbinom(n = length(prob_contam_fr_patdis), size = 1, prob = prob_contam_fr_patdis)
          # )
          ungroup() |>
          select(rid_uniq, prob_contam)

        ## step 9e: update room contam status #########################################
        ## room contam status updated based on 9c & 9d
        new_contam_stat = rbinom(
          n = nrow(room_contam_probs),
          size = 1,
          prob = room_contam_probs$prob_contam
        )
        idx_new_contam = which(new_contam_stat == 1)
        room_list$contam_next[idx_new_contam] = 1L
        # room_list$contam[idx_new_contam] = 1L

        ## step 9f: update HCW contam based on room contam ############################
        idx_room_contam = (room_list$contam@i + 1)
        hcw_contam_probs = tb_hcw_room_inter |>
          mutate(
            prob_hcw_now_contam = if_else(rid_uniq %in% idx_room_contam, (1 - prop_soap_out) * (1 - exp(-total_sec * tau_room_hcw)), 0)
          ) |>
          select(hid_uniq, prob_hcw_now_contam) |>
          filter(prob_hcw_now_contam > 0) |>
          group_by(hid_uniq) |>
          summarise(
            prob_hcw_not_contam_anyroom = prod(1 - prob_hcw_now_contam, na.rm = TRUE),
            .groups = "drop"
          ) |>
          mutate(prob_c = 1 - prob_hcw_not_contam_anyroom)

        new_contam_stat = rbinom(n = nrow(hcw_contam_probs), size = 1, prob = hcw_contam_probs$prob_c)
        idx_new_contam = hcw_contam_probs$hid_uniq[which(new_contam_stat == 1)]
        hcw_list$contam_next[idx_new_contam] = 1L

        ## step 9g: update HCW disease stat ############################################
        ## HCW disease status updated based on time spent in contam room
          # idx_room_contam # from 9f
        hcw_asymp_prob = tb_hcw_room_inter |>
          mutate(
            prob_hcw_asymp_rm = if_else(rid_uniq %in% idx_room_contam, (1 - exp(-total_sec * tau_room_hcw)), 0)
          ) |>
          select(hid_uniq, prob_hcw_asymp_rm) |>
          filter(prob_hcw_asymp_rm > 0) |>
          group_by(hid_uniq) |>
          summarise(
            prob_hcw_not_asymp = prod(1 - prob_hcw_asymp_rm, na.rm = TRUE),
            .groups = "drop"
          ) |>
          mutate(prob_hcw_asymp = 1 - prob_hcw_not_asymp)

        new_dis_stat = rbinom(n = nrow(hcw_asymp_prob), size = 1, prob = hcw_asymp_prob$prob_hcw_asymp)
        idx_new_hcw_dis = hcw_asymp_prob$hid_uniq[which(new_dis_stat == 1)]
        hcw_list$colonized[idx_new_hcw_dis] = 1L

        ## step 9h: update patient disease stat ########################################
        ## patient disease status updated based on time in contam room
        pat_idx_in_contam = room_list$occup[idx_room_contam]
        pat_idx_in_contam = pat_idx_in_contam[which(pat_idx_in_contam != 0)]
        prob_pat_latent = (1 - exp(
          -1 * 21600 * ## 21,600 sec in 6 hrs
            age_lambdas[pat_list$age_cat[pat_idx_in_contam]] * tau_room_patge65
        ))
        new_latent_stat = rbinom(
          n = length(prob_pat_latent),
          size = 1,
          prob = prob_pat_latent
        )
        idx_new_pat_latent_temp = pat_idx_in_contam[which(new_latent_stat == 1)]
        ## make sure "new" latent patients are not in other disease states, if so remove
        length(idx_new_pat_latent_temp)
        other_states = c((latent@i + 1),
                         (incubation@i + 1),
                         (asymptomatic@i + 1),
                         (symptomatic@i + 1),
                         (clin_res@i + 1),
                         (death@i + 1)
        )
        which(idx_new_pat_latent_temp %in% other_states) |> length()
        which(!idx_new_pat_latent_temp %in% other_states) |> length()
        idx_new_pat_latent_temp[!(idx_new_pat_latent_temp %in% other_states)] |> length()
        idx_new_pat_latent_temp[which(!idx_new_pat_latent_temp %in% other_states)] |> length()
        ## make sure "new" latent patients are in only coming from "susceptible" disease states, if not then remove
        idx_new_pat_latent_temp[idx_new_pat_latent_temp %in% (latent@i + 1)] |> length()
        idx_new_pat_latent = idx_new_pat_latent_temp[idx_new_pat_latent_temp %in% (susceptible@i + 1)]

        ## set number of days in incubation period for newly infected patients
        latent[idx_new_pat_latent] = 2L
        ## new latent period patients are no longer susceptible
        susceptible[idx_new_pat_latent] = 0L
        ### NEED to update/alter this portion, with tau parameter

        # ## 50/50 if incubation period is 1 or 2 days
        # rbinom(n = 10, size = 1, prob = c(0.5, 0.5)) + 1
        # table(rbinom(n = 100, size = 1, prob = c(0.5, 0.5)) + 1)

        ## reset contam status for rooms
        idx_room_contam = (room_list$contam_next@i + 1)
        room_list$contam[] = 0L
        room_list$contam[idx_room_contam] = 1L
        room_list$contam_next[] = 0L
        ## reset contam status for hcws
        idx_hcw_contam = (hcw_list$contam_next@i + 1)
        hcw_list$contam[] = 0L
        hcw_list$contam[idx_hcw_contam] = 1L
        hcw_list$contam_next[] = 0L
      } ## end of time block loop

      ## step 10: pat days left in stay ###############################################
      # for patients with a next visit & who are NOT in the hospital, decrease their # of days to next visit
      pat_w_next = pat_list$days_to_next_viz@i + 1
      pat_w_next_nihosp = pat_w_next[!pat_w_next %in% (room_list$occup@x)]
      idx_to_dec = which((pat_w_next) %in% pat_w_next_nihosp) # indices specific to the sparse matrix not patient index
      pat_list$days_to_next_viz@x[idx_to_dec] = (pat_list$days_to_next_viz@x[idx_to_dec] - 1)
      # ^ if a patient has a sympt. viz during their inter-viz time window, the in b/n time is essentially paused

      ## update patient days left in stay (moved from step 4)
      ## patients did all time blocks and are now ready for discharge
      ## indices where someone just completed last day
      idx_to_discharge = (pat_list$days_rem@i + 1)[which(pat_list$days_rem@x == 1)]
      ## decrease days rem for patients still in hosp & essentially discharge patients at this var
      pat_list$days_rem@x = pat_list$days_rem@x - 1

      ## update observed cdi case results
      symp_pat_idx = symptomatic@i + 1
      symp_in_hosp = symp_pat_idx[symp_pat_idx %in% (room_list$occup@x)]
      new_obs_symp = sum(!symp_in_hosp %in% symp_in_hosp_yesterday)
      symp_in_hosp_yesterday = symp_in_hosp          # symp_in_hosp_yesterday defined in initialization
      obs_results$obs_all_symp = obs_results$obs_all_symp + new_obs_symp
      ## update observed recurrent case resutls
      # recur_pat_idx & recur_pat_today  # from earlier in code
      recur_in_hosp = recur_pat_today[recur_pat_today %in% (room_list$occup@x)]
      new_obs_recur = sum(!recur_in_hosp %in% recur_in_hosp_yesterday)
      recur_in_hosp_yesterday = recur_in_hosp
      obs_results$obs_recur = obs_results$obs_recur + new_obs_recur

      tot_pat_newsymp_day[d] = obs_results$obs_all_symp

      ## calculate prevalence
      if(d == 32) {
        obs_results$prev_feb = length(symptomatic@i) / length((pat_list$days_rem)@i)
      }
      if(d == 60) {
        obs_results$prev_mar = length(symptomatic@i) / length((pat_list$days_rem)@i)
      }
      if(d == total_days_sim) {
        obs_results$prop_rm_contam_end = length(room_list$contam@x) / length(room_list$rid_uniq)
      }

    } ## end of day loop
    end_time = Sys.time()
    time_taken = end_time - start_time
    print(time_taken)
    to_return = list(
      tot_pat_newsymp_day = tot_pat_newsymp_day,
      prev_feb = obs_results$prev_feb,
      prev_mar = obs_results$prev_mar,
      prop_rm_contam_end = obs_results$prop_rm_contam_end
    )
    return(to_return)
    # return(tot_pat_newsymp_day)
  }
