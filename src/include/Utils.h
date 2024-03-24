#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include <array>

std::vector<std::string> split(const std::string& s, char delimiter);  //splitting string into chunks by given charachter
float calculate_distance(std::array<float,3> point_a, std::array<float,3> point_b); // calculates distance between two given points
std::vector <float> calculate_vector(float ax, float ay, float az, float bx, float by, float bz);
float calculate_dot_product(std::vector<float> vect_a, std::vector<float> vect_b);
std::vector<float> calculate_cross_product(std::vector<float> vect_a, std::vector<float> vect_b);

#endif