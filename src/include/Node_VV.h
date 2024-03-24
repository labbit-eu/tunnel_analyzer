#ifndef Node_VV_H
#define Node_VV_H

#include "Point_3d.h"
#include "Node_QTV.h"
#include "Gate.h"
#include <vector>

class Node_VV : public Point_3d{
    private:
        int vv_id_;
        std::vector<Node_QTV*> tangent_qtvs_;
        Gate all_gates_[4];
        std::vector<Node_VV*> connected_vv_nodes_;
        std::vector<float> distances_to_connected_nodes_;
        bool is_discarded_;
        bool is_starting_point_;
        bool is_user_starting_point_;
        bool is_alpha_starting_point_;
        bool is_minima_;
        bool is_visited_;
        bool is_in_priority_queue_;
        std::string dbscan_label_;
        float tentative_distance_;

    public:
        Node_VV();
        Node_VV(int n_id, float n_x, float n_y, float n_z);
        Node_VV(int n_id, float n_x, float n_y, float n_z, float n_radius);
        Node_VV(int n_id, float n_x, float n_y, float n_z, float n_radius, bool minima);

        void set_tangent_qtv_at(Node_QTV* t_qtv);
        void update_radius();
        void create_gates();
        void add_connection(Node_VV* node);
        void swap_connection(int node_id_to_swap, Node_VV* node);
        void remove_connection(int id_to_remove);
        void calculate_distance_to_neighbours();

        int get_id();
        Gate* get_gates();
        std::vector<Node_VV*> get_connections();
        std::vector<float> get_distances_to_neighbours();
        bool get_is_discarded();
        bool get_is_starting_point();
        bool get_is_user_starting_point();
        float get_tentative_distance();
        bool get_is_visited();
        bool get_is_in_priority_queue();
        bool get_is_minima();
        bool get_is_alpha_starting_point();
        std::string get_dbscan_label();
        std::vector<Node_QTV*> get_tangent_qtvs();

        void set_dbscan_label(std::string label);
        void set_is_discarded(bool arg);
        void set_is_starting_point(bool arg);
        void set_is_user_starting_point(bool arg);
        void set_tentative_distance(float value);
        void set_is_visited(bool arg);
        void set_is_in_priority_queue(bool arg);
        void set_is_minima(bool arg);
        void set_is_alpha_starting_point(bool arg);

        std::string print_to_string();
        std::string print_to_pdb_string();
};

#endif