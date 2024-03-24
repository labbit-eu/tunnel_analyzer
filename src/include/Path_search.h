#ifndef PATH_SEARCH
#define PATH_SEARCH

#include "Voronoi_diagram.h"
#include <vector>
#include <map>

class Path_search{
    private:
        std::vector<Node_VV> *voronoi_diagram_vertices_;
        std::vector<Node_VV*> starting_points_;
        std::vector<Node_VV*> user_starting_points_;
        std::map<int, Node_VV> indexed_nodes_;

    public:
        Path_search(std::vector<Node_VV> *vertices);

        void save_starting_points_to_vector();
        void closest_common_starting_point(std::array<float, 3> qt_center);
        void discard_starting_point_with_no_connections();
        int clusterize_starting_points(float epsilon, int minimal_no_points_in_cluster);
        void save_starting_points_to_pdb(std::string file_name);
        void save_starting_points_to_pdb(std::string file_name, std::vector<Node_VV *> points_vector);
        void save_shortest_path_tree_to_pdb(std::map<int, int> results, std::vector<Node_VV> *voronoi, std::string location = "");

        std::map<int,int> simple_dijkstra_path_search(std::vector<Node_VV> *voronoi, Node_VV* initial);
        void multiple_points_path_search();
        void multiple_points_path_search(std::vector<Node_VV *> user_stated_starting_points);

        std::vector<Node_VV*> get_starting_points();
        std::vector<Node_VV*> get_user_starting_points();
};

#endif