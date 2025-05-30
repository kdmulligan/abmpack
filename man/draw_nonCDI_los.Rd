% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_nonCDI_los.R
\name{draw_nonCDI_los}
\alias{draw_nonCDI_los}
\title{draw non-CDI los before simulation}
\usage{
draw_nonCDI_los(dist = noCDI_los_dist_mdc_tran6days, transfer, md_cat)
}
\arguments{
\item{dist}{a list object with los distributions by mdc and transfer type}

\item{transfer}{transfer type of the hospitalization: not a transfer,
transfer not last, transfer last}

\item{md_cat}{major diagnostic category, 0-25}
}
\value{
Returns a integer (or vector, function is vectorized) with the new LOS
}
\description{
If a patient has CDI according to the HCUP data, their LOS
needs to be redrawn from the non_CDI LOS. The new LOS is drawn from the
matching transfer/MDC distribution with days between 1-365
}
