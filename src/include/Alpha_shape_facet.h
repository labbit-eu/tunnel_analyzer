#ifndef ALPHA_SHAPE_FACET
#define ALPHA_SHAPE_FACET

#include <array>
#include <vector>

class Alpha_shape_facet{
    private:
        int face_id_;
        std::array<std::array<float, 3>, 3> vertices_of_face_;
        std::array<std::array<float,3>,2> first_edge_;
        std::array<std::array<float,3>,2> second_edge_;
        std::array<std::array<float,3>,2> third_edge_;
        std::vector<float> normal_;
        float constant_d_;

    public:
        Alpha_shape_facet(int id, std::array<std::array<float, 3>, 3> vertices);

        void calculate_d();
        void calculate_normal();
        bool check_if_neighbours(Alpha_shape_facet target_facet);
        std::array<std::array<float, 3>, 3> get_vertices_of_face();
        std::array<std::array<std::array<float,3>,2>,3> get_edges();
        std::vector<float> get_normal();
        int get_id();
        std::string print_to_pdb_string();
};

#endif