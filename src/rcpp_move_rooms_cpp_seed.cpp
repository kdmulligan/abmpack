//' @name move_rooms_cpp_seed
//' @title Move patient rooms
//' @description Move patients from one room to another (similar to assign
//' rooms function). Designed to deal with one facility at a time.
//' @param pat_rm_type Data frame of patients and the room type they are to go to.
//' @param icu vector of number of icu rooms available at each hospital
//' @param non vector of number of non-icu rooms available at each hospital
//' @param seed seed to be passed in to rcpp
//' @return returns data frame with one column for patient and a second column
//' for the room number to which they are assigned.

#include <Rcpp.h>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <random>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
Rcpp::DataFrame move_rooms_cpp_seed(DataFrame pat_rm_type, SEXP icu, SEXP non, unsigned int seed) {
  // `pat` df columns are visitlink, ahaid, adrgriskmortality
  int n = pat_rm_type.nrow();
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

  NumericVector pat_id = pat_rm_type[0];
  NumericVector room_type = pat_rm_type[1];

  // Initialize random number generator with seed
  std::mt19937 rng(seed);

  // Create a vector to store assigned rooms
  Rcpp::IntegerVector assignedRooms(pat_id.size());

  for(int i = 0; i < n; i++) {
    if(room_type[i] == 1) {   //room type is icu
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
        // Randomly select an ICU available room
        std::uniform_int_distribution<std::size_t> dist(0, non_rooms.size() - 1);
        std::size_t roomIndex = dist(rng);
        // Assign the room to the patient
        assignedRooms[i] = non_rooms[roomIndex];
        // Remove the assigned room from the vector of available rooms
        non_rooms.erase(non_rooms.begin() + roomIndex);
      }
    } else {  //room type is non-icu
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
