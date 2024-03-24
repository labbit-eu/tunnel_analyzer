#ifndef QUASI_TRIANGULATION
#define QUASI_TRIANGULATION

#include "Node_QTV.h"
#include "Node_QTC.h"
#include "Alpha_shape_facet.h"
#include <vector>
#include <list>
#include <fstream>
#include <array>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;
typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;
typedef Gt::Point_3                                  Point;
typedef Alpha_shape_3::Alpha_iterator                Alpha_iterator;

// class containing qt vertices and cells

class Quasi_triangulation{
    friend class Instance;
    private:
        std::list<Point> cgal_quasi_triangulation_points_; // necessary for CGAL 3d-alpha_shape calculations
        std::vector<Node_QTV> quasi_triangulation_vertices_; // qt vertices from .qtf file
        std::vector<Node_QTC> quasi_triangulation_cells_; // qt cells from .qtf file
        std::vector<Alpha_shape_facet> concave_hull_facets_; // concave hull facets acquired from cgal 3d alpha shape calculation
        std::vector<Alpha_shape_facet> larger_concave_hull_facets_;
        int alpha_value_;
        std::array<float, 3> center_coordinates_of_qt_; //average of  x, y, z coordinates
        std::array<std::pair<float, float>, 3> extreme_coordinates_; // extreme coordinates of protein's atoms [[x, -x],[y, -y],[z, -z]]

    public:
        void add_cgal_point(float x, float y, float z);
        void add_vertex(Node_QTV vertex);
        void add_cell(Node_QTC cell);
        void initialize_center_coordinates(std::array<float, 3> xyz_sums);
        void set_extreme_coordinates(float x, float neg_x, float y, float neg_y, float z, float neg_z);

        void calculate_alpha_shape();

        std::vector<Node_QTV> get_qt_vertices();
        std::vector<Node_QTC> get_qt_cells();
        std::array<float, 3> get_qt_center();
        std::vector<Alpha_shape_facet> get_concave_hull();
        std::vector<Alpha_shape_facet> get_larger_concave_hull();
        std::array<std::pair<float, float>, 3> get_extreme_coordinates();
        float calculate_extreme_radius();

        void save_alpha_shape_to_pdb(std::string file_name, Alpha_shape_3* alpha_shape);
        void save_my_facets_to_pdb();
};

#endif