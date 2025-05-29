#include <Rcpp.h>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <random>

using namespace Rcpp;
using namespace std;

//' @name rcpp_assign_rooms_cpp_seed
//' @title Assign rooms to patients
//' @description Randomly assign patients to hospital rooms with rcpp seed.
//' Designed to deal with one facility at a time.
//' @param pat_risks data frame where each row is a facility and four columns
//' one for each risk level to hold the number of patients at each risk level
//' per hospital
//' @param icu vector of number of icu rooms available at each hospital
//' @param non vector of number of non-icu rooms available at each hospital
//' @param seed seed to be passed in to rcpp
//' @return returns data frame with one column for patient and a second column
//' for the room number to which they are assigned.
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_assign_rooms_cpp_seed(DataFrame pat_risks, SEXP icu, SEXP non, unsigned int seed) {

  // `pat` df columns are visitlink, ahaid, adrgriskmortality
  int n = pat_risks.nrow();
  IntegerVector icu_rooms;
  IntegerVector non_rooms;
  // case to handle if facility has no ICU rooms at all or any available
  if (Rf_isNull(icu)) {
    icu_rooms = IntegerVector(); // empty integer vector if input is NULL
  } else {
    icu_rooms = as<IntegerVector>(icu); // Convert input to IntegerVector
  }
  // case to handle if facility has no NON rooms available
  if (Rf_isNull(non)) {
    non_rooms = IntegerVector(); // empty integer vector if input is NULL
  } else {
    non_rooms = as<IntegerVector>(non); // Convert input to IntegerVector
  }

  NumericVector pat_id = pat_risks[0];
  NumericVector mort_risk = pat_risks[1];

  // Initialize random number generator
  std::mt19937 rng(seed);

  // Create a vector to store assigned rooms
  Rcpp::IntegerVector assignedRooms(pat_id.size());

  for(int i = 0; i < n; i++) {
    if(mort_risk[i] == 4) {
      if(icu_rooms.size() > 0) {
        // Randomly select an ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, icu_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = icu_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        icu_rooms.erase(icu_rooms.begin() + roomIndex);
      } else{
        // draw from non vector
        // Randomly select an NON-ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, non_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = non_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        non_rooms.erase(non_rooms.begin() + roomIndex);
      }
    }
    if(mort_risk[i] == 3) {
      if(non_rooms.size() > 0) {
        // draw from non vector
        // Randomly select a NON-ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, non_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = non_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        non_rooms.erase(non_rooms.begin() + roomIndex);
      } else{
        // Randomly select an ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, icu_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = icu_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        icu_rooms.erase(icu_rooms.begin() + roomIndex);
      }
    }
    if(mort_risk[i] == 2) {
      if(non_rooms.size() > 0) {
        // draw from non vector
        // Randomly select a NON-ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, non_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = non_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        non_rooms.erase(non_rooms.begin() + roomIndex);
      } else{
        // Randomly select an ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, icu_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = icu_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        icu_rooms.erase(icu_rooms.begin() + roomIndex);
      }
    }
    if(mort_risk[i] == 1) {
      if(non_rooms.size() > 0) {
        // draw from non vector
        // Randomly select a NON-ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, non_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = non_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        non_rooms.erase(non_rooms.begin() + roomIndex);
      } else{
        // Randomly select an ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, icu_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = icu_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        icu_rooms.erase(icu_rooms.begin() + roomIndex);
      }
    }
  }
  // Create a new data frame with the assigned rooms
  return Rcpp::DataFrame::create(Rcpp::Named("patid") = pat_id,
                                 Rcpp::Named("assigned_room") = assignedRooms);

}
