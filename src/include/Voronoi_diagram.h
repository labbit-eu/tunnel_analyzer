#ifndef VORONOI_DIAGRAM
#define VORONOI_DIAGRAM

#include "Node_VV.h"
#include "Alpha_shape_facet.h"
#include "Utils.h"
#include <vector>


//class containing voronoi diagram

class Voronoi_diagram{
    friend class Instance;
    protected:
        std::vector<Node_VV> voronoi_diagram_vertices_;
        int voronoi_diagram_initial_size_;

    public:
        void add_vertex(Node_VV vertex);
        void trim_discarded_vertices();
        void filter_by_user_radius(float user_radius_filter_value);
        void filter_vertices_outside_of_concave_hull(std::vector<Alpha_shape_facet> concave_hull, std::array<float, 3> center);
        void calculate_starting_points_from_alpha_shapes(std::vector<Alpha_shape_facet> concave_hull, std::vector<Alpha_shape_facet> larger_concave_hull);
        void calculate_distances_to_neighbours();
        void calculate_minimas(int user_radius_filter_value);
        void mark_user_starting_points(float search_radius, std::vector<std::array<float,3>> user_points);
        void search_for_user_starting_points(float search_radius, std::vector<std::array<float,3>> user_points);

        void set_initial_voronoi_size();
        void save_minimas_to_pdb(std::string file_name);
        void save_diagram_to_pdb(std::string file_name, bool save_all = true);
        void save_diagram_to_mmcif_with_radius(std::string file_name, bool save_all = true);
        void save_diagram_to_mmcif(std::string file_name, bool save_all = true);
        void save_diagram_to_output_json();
        void save_starting_points_to_pdb(std::string file_name);
        std::vector<Node_VV> get_voronoi_diagram();
        std::vector<Node_VV> *get_voronoi_diagram_adress();
};

#endif