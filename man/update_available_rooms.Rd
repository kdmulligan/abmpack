% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_available_rooms.R
\name{update_available_rooms}
\alias{update_available_rooms}
\title{Update list of available rooms based on occupancy}
\usage{
update_available_rooms(rm_list, hcupids_vec)
}
\arguments{
\item{rm_list}{list object with room info}

\item{hcupids_vec}{vector of hcup ids}
}
\value{
Data frame with number of rows equal to the number of patient moves per day
}
\description{
This function has no parameters because it uses the room_list
object from the environment. Additionally, it has no return object because
it replaces icu_avail and non_avail in the environment.
}
