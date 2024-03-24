#include "../include/Alpha_shape_facet.h"
#include "../include/Utils.h"
#include <iomanip>

Alpha_shape_facet::Alpha_shape_facet(int id, std::array<std::array<float, 3>, 3> vertices){
    face_id_ = id;
    vertices_of_face_ = vertices;
    calculate_normal();
    calculate_d();
    first_edge_ = {vertices[0], vertices[1]};
    second_edge_ = {vertices[0], vertices[2]};
    third_edge_ = {vertices[1], vertices[2]};
}

void Alpha_shape_facet::calculate_d(){
	constant_d_ = -(normal_[0] * vertices_of_face_[0][0]) - (normal_[1] * vertices_of_face_[0][1]) - (normal_[2] * vertices_of_face_[0][2]);
}

void Alpha_shape_facet::calculate_normal(){
	normal_ = calculate_cross_product(calculate_vector(vertices_of_face_[0][0], vertices_of_face_[0][1], vertices_of_face_[0][2], vertices_of_face_[1][0], vertices_of_face_[1][1], vertices_of_face_[1][2]),
                                     calculate_vector(vertices_of_face_[0][0], vertices_of_face_[0][1], vertices_of_face_[0][2], vertices_of_face_[2][0], vertices_of_face_[2][1], vertices_of_face_[2][2]));
}

bool Alpha_shape_facet::check_if_neighbours(Alpha_shape_facet target_facet){
    if(this->get_id() == target_facet.get_id()){
        return false;
    }else{
        int matches_count = 0;
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                if(this->get_vertices_of_face()[i] == target_facet.get_vertices_of_face()[j]){
                    matches_count ++;
                }
            }
        }
        if(matches_count == 2){
            return true;
        }else{
            return false;
        }
    }
}

std::array<std::array<float, 3>, 3> Alpha_shape_facet::get_vertices_of_face(){
    return vertices_of_face_;
}

std::vector<float> Alpha_shape_facet::get_normal(){
	return normal_;
}

int Alpha_shape_facet::get_id(){
    return face_id_;
}

std::array<std::array<std::array<float,3>,2>,3> Alpha_shape_facet::get_edges(){
    return {first_edge_, second_edge_, third_edge_};
}

std::string Alpha_shape_facet::print_to_pdb_string(){
    std::ostringstream output;
    output          << "ATOM  " << std::setw(5) << 5
                        << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[0][0] << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[0][1] << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[0][2] << "      " << std::setw(6) << std::setprecision(3)
                        << 1 << "           H  " << std::endl;
    output          << "ATOM  " << std::setw(5) << 5
                        << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[1][0] << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[1][1] << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[1][2] << "      " << std::setw(6) << std::setprecision(3)
                        << 1 << "           H  " << std::endl;
    output          << "ATOM  " << std::setw(5) << 5
                        << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[2][0] << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[2][1] << std::setw(8) << std::setprecision(3)
                        << this->get_vertices_of_face()[2][2] << "      " << std::setw(6) << std::setprecision(3)
                        << 1 << "           H  " << std::endl;
    return output.str();
}