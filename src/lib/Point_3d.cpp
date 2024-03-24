#include "../include/Point_3d.h"

Point_3d::Point_3d(){
    x_ = 999.999;
    y_ = 999.999;
    z_ = 999.999;
    radius_ = 999.999;
}

Point_3d::Point_3d(float x, float y, float z){
    x_ = x;
    y_ = y;
    z_ = z;
    radius_ = 999.999;
}

Point_3d::Point_3d(float x, float y, float z, float radius){
    x_ = x;
    y_ = y;
    z_ = z;
    radius_ = radius;
}

float Point_3d::get_x(){
    return x_;
}

float Point_3d::get_y(){
    return y_;
}

float Point_3d::get_z(){
    return z_;
}

float Point_3d::get_radius(){
    return radius_;
}

std::array<float,3> Point_3d::get_coordinates(){
    std::array<float,3> coordinates_array = {x_, y_, z_};
    return coordinates_array;
}

void Point_3d::set_x(float x){
    x_ = x;
}

void Point_3d::set_y(float y){
    y_ = y;
}

void Point_3d::set_z(float z){
    z_ = z;
}

void Point_3d::set_radius(float radius){
    radius_ = radius;
}
