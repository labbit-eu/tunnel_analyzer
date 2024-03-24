#include "../include/Utils.h"


std::vector<std::string> split(const std::string& s, char delimiter){
    // function for splitting string by given sign
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)){
        tokens.push_back(token);
    }
    return tokens;
}



float calculate_distance(std::array<float,3> point_a, std::array<float,3> point_b){
    // calculates distance between two given points
    float x1 = pow((point_a[0] - point_b[0]), 2);
    float y1 = pow((point_a[1] - point_b[1]), 2);
    float z1 = pow((point_a[2] - point_b[2]), 2);
    return sqrt(x1 + y1 + z1);
}

std::vector <float> calculate_vector(float ax, float ay, float az, float bx, float by, float bz){
    /*  function calculating a vector between two points
        input: 6 floats -> coordinates x,y,z of point a and b
        output: return vector structure -> calculated vector
    */
    std::vector<float> my_vector;
    my_vector.push_back(bx - ax);
    my_vector.push_back(by - ay);
    my_vector.push_back(bz - az);
    return my_vector;
}

float calculate_dot_product(std::vector<float> vect_a, std::vector<float> vect_b){
    float dot = ((vect_a[0]*vect_b[0]) + (vect_a[1]*vect_b[1]) + (vect_a[2]*vect_b[2]));
    return dot;
}

std::vector<float> calculate_cross_product(std::vector<float> vect_a, std::vector<float> vect_b){
    std::vector<float> cross;
    cross.push_back((vect_a[1]*vect_b[2]) - (vect_a[2]*vect_b[1]));
    cross.push_back((vect_a[2]*vect_b[0]) - (vect_a[0]*vect_b[2]));
    cross.push_back((vect_a[0]*vect_b[1]) - (vect_a[1]*vect_b[0]));
    return cross;
}
