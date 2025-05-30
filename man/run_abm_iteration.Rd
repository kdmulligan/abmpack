% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_abm_iteration.R
\name{run_abm_iteration}
\alias{run_abm_iteration}
\title{Run one ABM simulation}
\usage{
run_abm_iteration(
  n_days = 72,
  t_symp_room = 5e-04,
  t_hcw_room = 5e-04,
  t_room_hcw = 5e-04,
  t_room_patge65 = 5e-04,
  t_prob_room_contam_init = 0.005,
  l_poi_latent = 1,
  ind_col_pat_trans = 1,
  ind_asymp_pat_trans = 1,
  ind_asymp_hcw_trans = 1,
  ind_recur_trans = 1,
  ind_transfer_pat_trans = 1,
  SEED = 1212
)
}
\arguments{
\item{n_days}{number of days for the simulation}

\item{t_symp_room}{transmission from symp pat to room prior value}

\item{t_hcw_room}{transmission from hcw to room prior value}

\item{t_room_hcw}{transmission from room to hcw prior value}

\item{t_room_patge65}{transmission from room to patient 65 or older prior value}

\item{t_prob_room_contam_init}{probability a room is contam prior value}

\item{l_poi_latent}{poisson lambda value for latent per prior value}

\item{ind_col_pat_trans}{indicator for colonized patient transmission}

\item{ind_asymp_pat_trans}{indicator for asymptomatic patient transmission}

\item{ind_asymp_hcw_trans}{indicator for asymptomatic hcw transmission}

\item{ind_recur_trans}{indicator for recurrent patient transmission}

\item{ind_transfer_pat_trans}{indicator for transfer patient transmission}

\item{SEED}{seed for the simulation}
}
\value{
numeric value with total number of observed cases
}
\description{
Take a input of parameters and indicators to run simulation by
}
