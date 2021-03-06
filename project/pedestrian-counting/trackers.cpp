/* trackers.cpp
 * a class that defines tracking functions for the pedestrian couter
 * note that all methods in this class will be static
 
 * Brian J Gravelle
 * ix.cs.uoregon.edu/~gravelle
 * gravelle@cs.uoregon.edu

 * If by some miracle you find this software useful, thanks are accepted in
 * the form of chocolate or introductions to potential employers.

 */

#include "trackers.h"

double Trackers::max_dist_sqd = 1000;
char   Trackers::algo         = 'D';

Object* Trackers::find_previous_object(vector<Object> &old_objs, Object &curr_obj) {
  switch(algo) {
    case 'O':
      return find_previous_object_overlap(old_objs, curr_obj);
      break;
    case 'D':
    default:
      return find_previous_object_dist(old_objs, curr_obj);
  }
  return NULL;
}

Object* Trackers::find_previous_object_dist(vector<Object> &old_objs, Object &curr_obj) {
  double dist, min_dist = -1.0;
  Object *prev_obj = NULL;
  for(vector<Object>::iterator it_old_obj = old_objs.begin(); it_old_obj != old_objs.end(); it_old_obj++) {
    dist = it_old_obj->find_distance_sqd(curr_obj);
    if( (dist <= max_dist_sqd) && ((min_dist < 0) || (dist < min_dist)) ) {
      min_dist = dist;
      prev_obj = &(*it_old_obj); 
    }
  } 
  return prev_obj;
}

Object* Trackers::find_previous_object_overlap(vector<Object> &old_objs, Object &curr_obj) {
  double area, max_area = 0.0;
  Object *prev_obj = NULL;
  for(vector<Object>::iterator it_old_obj = old_objs.begin(); it_old_obj != old_objs.end(); it_old_obj++) {
    area = it_old_obj->find_overlap_area(curr_obj);
    if( area > max_area ) {
      max_area = area;
      prev_obj = &(*it_old_obj);
    }
  }
  return prev_obj;
}


char Trackers::get_algo(){
  return algo;
}

void Trackers::set_algo(char new_algo){
  algo = new_algo;
}

double Trackers::get_max_dist_sqd(){
  return max_dist_sqd;
}

void   Trackers::set_max_dist_sqd(double new_max_dist_sqd){
  max_dist_sqd = new_max_dist_sqd;
}


