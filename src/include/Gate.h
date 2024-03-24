#ifndef GATE_H
#define GATE_H

#include "Node_QTV.h"
#include <vector>

class Node_VV;

class Gate{
    private:
        int gate_id_;
        float normal_vector_[3], constant_d_, barycenter_[3];
        std::vector<Node_QTV*> gate_qtv_nodes_;
        Node_QTV* cell_qtv_node_;
        Node_VV* start_vv_node_;
        Node_VV* end_vv_node_;

    public:
        Gate();
        Gate(int n_id, Node_VV* s_vv_node, Node_QTV* n_qtv1, Node_QTV* n_qtv2, Node_QTV* n_qtv3);
        Gate(int n_id, Node_VV* s_vv_node, Node_QTV* n_qtv1, Node_QTV* n_qtv2, Node_QTV* n_qtv3, Node_QTV* n_cell);

        void set_end_vv_node(Node_VV* end_node);
        Node_VV* get_start_vv_node();
        Node_VV* get_end_vv_node();
        std::vector<Node_QTV*> get_gate_qtv_nodes();
        std::string print_to_string();
};

#endif