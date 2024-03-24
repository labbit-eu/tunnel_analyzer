#include "../include/Gate.h"
#include "../include/Node_VV.h"

#include <iomanip>
#include <iostream>

class Node_VV;
Gate::Gate(){
    gate_id_ = -1;
}

Gate::Gate(int n_id, Node_VV* s_vv_node, Node_QTV* n_qtv1, Node_QTV* n_qtv2, Node_QTV* n_qtv3){
    gate_id_ = n_id;
    start_vv_node_ = s_vv_node;
    end_vv_node_ = new Node_VV;
    gate_qtv_nodes_.push_back(n_qtv1);
    gate_qtv_nodes_.push_back(n_qtv2);
    gate_qtv_nodes_.push_back(n_qtv3);
}

Gate::Gate(int n_id, Node_VV* s_vv_node, Node_QTV* n_qtv1, Node_QTV* n_qtv2, Node_QTV* n_qtv3, Node_QTV* n_cell){
    gate_id_ = n_id;
    start_vv_node_ = s_vv_node;
    end_vv_node_ = new Node_VV;
    gate_qtv_nodes_.push_back(n_qtv1);
    gate_qtv_nodes_.push_back(n_qtv2);
    gate_qtv_nodes_.push_back(n_qtv3);
    cell_qtv_node_ = n_cell;
}

void Gate::set_end_vv_node(Node_VV *end_node){
    end_vv_node_ = end_node;
}

Node_VV* Gate::get_start_vv_node(){
    return start_vv_node_;
}

Node_VV* Gate::get_end_vv_node(){
    return end_vv_node_;
}

std::vector<Node_QTV*> Gate::get_gate_qtv_nodes(){
    return gate_qtv_nodes_;
}

std::string Gate::print_to_string(){
    std::cout << std::setprecision(3) << std::fixed;
    std::ostringstream output;
    output << "Gate: " << std::setw(2) << gate_id_ << std::endl;
    output << "\tId of start node: " << std::setw(3) << start_vv_node_->get_id() << std::endl;
    output << "\t\tStart Coordinates: " << std::setw(3) << start_vv_node_->get_x() << ", " << std::setw(3) << start_vv_node_->get_y() << ", "  << std::setw(3) << start_vv_node_->get_z() << std::endl;
    output << "\tIDs of QTV nodes: " << gate_qtv_nodes_[0]->get_qtv_atom_id() << ", "<< gate_qtv_nodes_[1]->get_qtv_atom_id() << ", "<< gate_qtv_nodes_[2]->get_qtv_atom_id() << std::endl;
    output << "\tID of the end node: " << std::setw(4) << end_vv_node_->get_id() << std::endl;
    output << "\t\tEnd Coordinates: " << std::setw(3) << end_vv_node_->get_x() << ", " << std::setw(3) << end_vv_node_->get_y() << ", "  << std::setw(3) << end_vv_node_->get_z() << std::endl;
    return output.str();
}