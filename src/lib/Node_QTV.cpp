#include "../include/Node_QTV.h"
#include <iomanip>
#include <iostream>

Node_QTV::Node_QTV():Point_3d(){
    qtv_id_ = -1;
    qtv_atom_id_ = -1;
}

Node_QTV::Node_QTV(int id, int atom_id, float x, float y, float z, float radius):Point_3d(x, y, z, radius){
    qtv_id_ = id;
    qtv_atom_id_ = atom_id;
}

int Node_QTV::get_qtv_id(){
    return qtv_id_;
}

int Node_QTV::get_qtv_atom_id(){
    return qtv_atom_id_;
}

std::string Node_QTV::print_to_pdb_string(){
    std::ostringstream output;
    output          << "ATOM  " << std::setw(5) << this->get_qtv_atom_id()
                    << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                    << this->get_x() << std::setw(8) << std::setprecision(3)
                    << this->get_y() << std::setw(8) << std::setprecision(3)
                    << this->get_z() << "      " << std::setw(6) << std::setprecision(3)
                    << this->get_radius() << "           H  " << std::endl;
    return output.str();
}