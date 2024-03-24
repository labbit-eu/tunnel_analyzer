#include "../include/Node_VV.h"
#include "../include/Utils.h"
#include <algorithm>
#include <iomanip>
#include <iostream>

Node_VV::Node_VV():Point_3d(){
    vv_id_ = -1;
    is_discarded_ = false;
    is_visited_ = false;
    is_starting_point_ = false;
    is_user_starting_point_ = false;
    tentative_distance_ = std::numeric_limits<float>::infinity();
    is_in_priority_queue_ = false;
    is_minima_ = false;
    is_alpha_starting_point_ = false;
    dbscan_label_ = "undefined";
}

Node_VV::Node_VV(int id, float x, float y, float z):Point_3d(x, y, z){
    vv_id_ = id;
    is_discarded_ = false;
    is_visited_ = false;
    is_starting_point_ = false;
    is_user_starting_point_ = false;
    tentative_distance_ = std::numeric_limits<float>::infinity();
    is_in_priority_queue_ = false;
    is_minima_ = false;
    is_alpha_starting_point_ = false;
    dbscan_label_ = "undefined";
}

Node_VV::Node_VV(int id, float x, float y, float z, float radius):Point_3d(x, y, z, radius){
    vv_id_ = id;
    is_discarded_ = false;
    is_visited_ = false;
    is_starting_point_ = false;
    is_user_starting_point_ = false;
    tentative_distance_ = std::numeric_limits<float>::infinity();
    is_in_priority_queue_ = false;
    is_minima_ = false;
    is_alpha_starting_point_ = false;
    dbscan_label_ = "undefined";
}

Node_VV::Node_VV(int id, float x, float y, float z, float radius, bool minima):Point_3d(x, y, z, radius){
    vv_id_ = id;
    is_discarded_ = false;
    is_visited_ = false;
    is_starting_point_ = false;
    is_user_starting_point_ = false;
    tentative_distance_ = std::numeric_limits<float>::infinity();
    is_in_priority_queue_ = false;
    is_minima_ = minima;
    is_alpha_starting_point_ = false;
    dbscan_label_ = "undefined";
}

void Node_VV::set_tangent_qtv_at(Node_QTV* t_qtv){
    tangent_qtvs_.push_back(t_qtv);
}

void Node_VV::update_radius(){
    float dist_a = calculate_distance(this->get_coordinates(), tangent_qtvs_[0]->get_coordinates());
    float dist_b = calculate_distance(this->get_coordinates(), tangent_qtvs_[1]->get_coordinates());
    float dist_c = calculate_distance(this->get_coordinates(), tangent_qtvs_[2]->get_coordinates());
    float dist_d = calculate_distance(this->get_coordinates(), tangent_qtvs_[3]->get_coordinates());
    float a_rad = dist_a - tangent_qtvs_[0]->get_radius();
    float b_rad = dist_b - tangent_qtvs_[1]->get_radius();
    float c_rad = dist_c - tangent_qtvs_[2]->get_radius();
    float d_rad = dist_d - tangent_qtvs_[3]->get_radius();
    this->set_radius((a_rad + b_rad + c_rad + d_rad) / 4.0);
}

void Node_VV::create_gates(){
    Gate gate1(1, this, tangent_qtvs_[1], tangent_qtvs_[2], tangent_qtvs_[3], tangent_qtvs_[0]);
    Gate gate2(2, this, tangent_qtvs_[0], tangent_qtvs_[2], tangent_qtvs_[3], tangent_qtvs_[1]);
    Gate gate3(3, this, tangent_qtvs_[0], tangent_qtvs_[1], tangent_qtvs_[3], tangent_qtvs_[2]);
    Gate gate4(4, this, tangent_qtvs_[0], tangent_qtvs_[1], tangent_qtvs_[2], tangent_qtvs_[3]);

    all_gates_[0] = gate1;
    all_gates_[1] = gate2;
    all_gates_[2] = gate3;
    all_gates_[3] = gate4;
}

void Node_VV::add_connection(Node_VV* node){
    bool exist_in = false;
    for (uint i = 0; i<connected_vv_nodes_.size(); i++){
        if (connected_vv_nodes_[i]->get_id() == node->get_id()){
            exist_in = true;
            break;
        }
    }
    if (!exist_in){
        connected_vv_nodes_.push_back(node);
    }
}

void Node_VV::swap_connection(int node_id_to_swap, Node_VV* node){
    for(int i = 0; i< connected_vv_nodes_.size();i++){
        if(connected_vv_nodes_[i]->get_id() == node_id_to_swap){
            connected_vv_nodes_[i] = node;
        }
    }
}

void Node_VV::remove_connection(int id_to_remove){
    connected_vv_nodes_.erase(remove_if(connected_vv_nodes_.begin(), connected_vv_nodes_.end(), [&id_to_remove](const Node_VV * s){ return s->vv_id_ == id_to_remove;}), connected_vv_nodes_.end());
}

void Node_VV::calculate_distance_to_neighbours(){
    for(uint i = 0; i < connected_vv_nodes_.size(); i++){
        float distance = calculate_distance(this->get_coordinates(), connected_vv_nodes_[i]->get_coordinates() );
        distances_to_connected_nodes_.push_back(distance);
    }
}

Gate* Node_VV::get_gates(){
    return all_gates_;
}

std::vector<Node_VV*> Node_VV::get_connections(){
    return connected_vv_nodes_;
}

std::vector<float> Node_VV::get_distances_to_neighbours(){
    return distances_to_connected_nodes_;
}

int Node_VV::get_id(){
    return vv_id_;
}

std::string Node_VV::get_dbscan_label(){
    return dbscan_label_;
}

void Node_VV::set_dbscan_label(std::string label){
    dbscan_label_ = label;
}

void Node_VV::set_is_discarded(bool arg){
    is_discarded_ = arg;
}

void Node_VV::set_is_user_starting_point(bool arg){
    is_user_starting_point_ = arg;
}

void Node_VV::set_is_starting_point(bool arg){
    is_starting_point_ = arg;
}

void Node_VV::set_tentative_distance(float value){
    tentative_distance_ = value;
}

void Node_VV::set_is_visited(bool arg){
    is_visited_ = arg;
}

void Node_VV::set_is_in_priority_queue(bool arg){
    is_in_priority_queue_ = arg;
}

void Node_VV::set_is_minima(bool arg){
    is_minima_ = arg;
}

void Node_VV::set_is_alpha_starting_point(bool arg){
    is_alpha_starting_point_ = arg;
}

bool Node_VV::get_is_minima(){
    return is_minima_;
}
bool Node_VV::get_is_discarded(){
    return is_discarded_;
}

bool Node_VV::get_is_user_starting_point(){
    return is_user_starting_point_;
}

bool Node_VV::get_is_starting_point(){
    return is_starting_point_;
}

float Node_VV::get_tentative_distance(){
    return tentative_distance_;
}

bool Node_VV::get_is_visited(){
    return is_visited_;
}

bool Node_VV::get_is_in_priority_queue(){
    return is_in_priority_queue_;
}

bool Node_VV::get_is_alpha_starting_point(){
    return is_alpha_starting_point_;
}

std::vector<Node_QTV*> Node_VV::get_tangent_qtvs(){
    return tangent_qtvs_;
}

std::string Node_VV::print_to_string(){
    std::ostringstream output;
    output  << "#------------------------------------------------------------------------------#" << std::endl;
    output  << "ID: " << vv_id_ << std::endl;
    output  << "\tPosition: " << std::setw(7) << this->get_x() << ", " << std::setw(7) << this->get_y() << ", "  << std::setw(7) << this->get_z() << std::endl;
    output  << "\tRadius: " << this->get_radius() << std::endl;
    output  << "\tIDs of tangent QTV_Nodes: " << std::setw(4) << tangent_qtvs_[0]->get_qtv_id() << ", " 
            << std::setw(4) << tangent_qtvs_[1]->get_qtv_id() << ", "
            << std::setw(4) << tangent_qtvs_[2]->get_qtv_id() << ", "
            << std::setw(4) << tangent_qtvs_[3]->get_qtv_id() << std::endl;

    output << "\tIDs of connected VV_Nodes: ";
    for (ushort i = 0; i < connected_vv_nodes_.size(); ++i){
        output << std::setw(4) << connected_vv_nodes_[i]->get_id() << ", ";
    }

    output << std::endl;
    output << "\tDistances to connected VV_Nodes: ";
    for(ushort i = 0; i < distances_to_connected_nodes_.size(); i++){
        output << std::setw(6) << distances_to_connected_nodes_[i] << ", ";
    }
    output << std::endl;
    output << "is discarded: "<< is_discarded_ <<std::endl;
    output << "is starting point: "<< is_starting_point_ <<std::endl;

    return output.str();
}

std::string Node_VV::print_to_pdb_string(){
    std::ostringstream output;
    output          << "ATOM  " << std::setw(5) << this->get_id()
                    << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                    << this->get_x() << std::setw(8) << std::setprecision(3)
                    << this->get_y() << std::setw(8) << std::setprecision(3)
                    << this->get_z() << "      " << std::setw(6) << std::setprecision(3)
                    << this->get_radius() << "           H  " << std::endl;
    return output.str();
}