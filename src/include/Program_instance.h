#ifndef PROGRAM_INSTANCE
#define PROGRAM_INSTANCE


#include "Voronoi_diagram.h"
#include "Quasi_triangulation.h"
#include "Path_search.h"

// single protein instance class

class Instance{
    private:
        std::string qtf_file_name_;
        float radius_filter_value_;
        Voronoi_diagram voronoi_diagram_;
        Quasi_triangulation quasi_triangulation_;
        std::vector<std::array<float,3>> user_points_of_search;

        void automatic_search();
        void single_starting_sphere_search();
        void multiple_starting_spheres_search();

    public:
        Instance(std::string snapshot_qtf_file_name);
        void parse_user_input_file(std::string user_input_file_name);
        void parse_qtf_file();
        void set_connections();
};

#endif