% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sample_day_mvts_cpp_seed}
\alias{sample_day_mvts_cpp_seed}
\title{Determine number of patient mvts for the day.}
\usage{
sample_day_mvts_cpp_seed(los, cur_room_type, seed)
}
\arguments{
\item{los}{Vector of patient los for the visit,}

\item{cur_room_type}{vector of patients' current room type (icu/non).}

\item{seed}{value of seed to be passed into the rcpp.}
}
\value{
returns character vector with the number of "icu,non" movements
for the patients' upcoming day.
}
\description{
function to randomly sample patient daily movement patterns
}
