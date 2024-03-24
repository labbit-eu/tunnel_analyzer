#include "../include/Path_search.h"
#include "../include/Utils.h"
#include <map>
#include <queue>
#include <limits>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;
extern std::ofstream result_json_file_alt;

Path_search::Path_search(std::vector<Node_VV> *voronoi_diagram){
    voronoi_diagram_vertices_ = voronoi_diagram;
}

void Path_search::save_starting_points_to_vector(){
    starting_points_.clear();
    for(int i = 0; i < voronoi_diagram_vertices_->size(); i++){
        if(!voronoi_diagram_vertices_->at(i).get_is_discarded()){
            indexed_nodes_[voronoi_diagram_vertices_->at(i).get_id()] = voronoi_diagram_vertices_->at(i);
            if(voronoi_diagram_vertices_->at(i).get_is_starting_point()){
                starting_points_.push_back(&voronoi_diagram_vertices_->at(i));
            }
            if(voronoi_diagram_vertices_->at(i).get_is_user_starting_point()){
                user_starting_points_.push_back(&voronoi_diagram_vertices_->at(i));
            }
        }
    }
}

struct CompareDistance{     //  struct used for sorting min priority queue
    bool operator()(Node_VV* v1, Node_VV* v2)
    {
        return v1->get_tentative_distance() > v2->get_tentative_distance();
    }
};

void Path_search::save_shortest_path_tree_to_pdb(std::map<int, int> results, std::vector<Node_VV> *voronoi, std::string location){
    /*  saves calculated shortest path tree in .pdb file for visualization purposes
    */
    std::ofstream file;
    file.open(location + "0_shortest_path_tree.pdb");
    for(uint i = 0; i < voronoi->size(); i++){
        if(!voronoi->at(i).get_is_discarded()){
            file << "ATOM  " << std::setw(5) << voronoi->at(i).get_id() << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                 << voronoi->at(i).get_x() << std::setw(8) << std::setprecision(3)
                 << voronoi->at(i).get_y() << std::setw(8) << std::setprecision(3)
                 << voronoi->at(i).get_z() << "      " << std::setw(6) << std::setprecision(3)
                 << voronoi->at(i).get_radius() << "           H  " << std::endl;
        }
    }
    for(std::map < int, int>::const_iterator it = results.begin(); it != results.end(); it++){
        if(!indexed_nodes_.at(it->first).get_is_discarded() && !indexed_nodes_.at(it->second).get_is_discarded()){
            file << "CONECT" << std::setw(5) << it->first << std::setw(5) << it->second << std::endl;
        }
    }
    file.close();
}

void print_total_cost(std::map<int, float> results, std::string log_file_name){
    std::ofstream file;
    file.open(log_file_name + ".txt");
    for(std::map<int, float>::const_iterator it = results.begin(); it != results.end(); it++){
        file << "Distance from initial node to " << it->first << " : " << it->second << std::endl;
    }
    file.close();
}

void print_previous_nodes(std::map<int, int> results, std::string log_file_name){
    std::ofstream file;
    file.open(log_file_name + ".txt");
    file << "Connections: " << std::endl;
    for(std::map<int, int>::const_iterator it = results.begin(); it != results.end(); it++){
        file << "from " << it->first << " to " << it->second << std::endl;
    }
    file.close();
}

std::map<int,int> Path_search::simple_dijkstra_path_search(std::vector<Node_VV> *voronoi, Node_VV* initial){
    std::map<int, float> total_cost;
    std::map<int, int> previous_nodes;
    std::priority_queue<Node_VV*, std::vector<Node_VV*>, CompareDistance> min_pq;

    for(uint i = 0; i<voronoi->size(); i++){
            voronoi->at(i).set_tentative_distance(std::numeric_limits<float>::infinity());
            voronoi->at(i).set_is_visited(false);
            total_cost[voronoi->at(i).get_id()] = voronoi->at(i).get_tentative_distance();
    }

    initial->set_tentative_distance(0); //start vertex has distance = 0, all other ones are already set to infinity
    initial->set_is_visited(false);
    total_cost[initial->get_id()] = initial->get_tentative_distance();

    min_pq.push(initial);
    bool path_found = false;
    int debug_test = 0;
    while(!min_pq.empty()){
        Node_VV *new_smallest = min_pq.top();
        min_pq.top()->set_is_visited(true);
        min_pq.top()->set_is_in_priority_queue(false);
        min_pq.pop();
        for(uint i = 0; i<new_smallest->get_connections().size(); i++){ //for each neighobur of closest node from pq
            if(!new_smallest->get_connections()[i]->get_is_discarded()){
                if(!new_smallest->get_connections()[i]->get_is_visited()){ //check if visited
                    if(!new_smallest->get_connections()[i]->get_is_in_priority_queue()){ //check if vertex is already in priority queue
                        //before adding neighbour to queue chek if it's already there, if it is, check the distance, if its distance is shorter, save it,
                        min_pq.push(new_smallest->get_connections()[i]); //push its neighbours to min_pq
                        new_smallest->get_connections()[i]->set_is_in_priority_queue(true);
                        float new_path_length = total_cost.find(new_smallest->get_id())->second + new_smallest->get_distances_to_neighbours()[i]; //calculate length from start node to this node
                        if(new_path_length < total_cost.find(new_smallest->get_connections()[i]->get_id())->second){ //if length is shorter
                            total_cost[new_smallest->get_connections()[i]->get_id()] = new_path_length; //update path length for neighbour
                            previous_nodes[new_smallest->get_connections()[i]->get_id()] = new_smallest->get_id(); //update previous node to node we removed
                            new_smallest->get_connections()[i]->set_tentative_distance(new_path_length); //update new shorter distance in min_pq
                            std::make_heap(const_cast<Node_VV**>(&min_pq.top()), const_cast<Node_VV**>(&min_pq.top())+ min_pq.size(), CompareDistance()); //sorts heap so that the closest by distance vertex is on top
                        }
                    }else{
                        float new_path_length = total_cost.find(new_smallest->get_id())->second + new_smallest->get_distances_to_neighbours()[i];
                        if(new_path_length < total_cost.find(new_smallest->get_connections()[i]->get_id())->second){ //if length is shorter
                            total_cost[new_smallest->get_connections()[i]->get_id()] = new_path_length; //update path length for neighbour
                            previous_nodes[new_smallest->get_connections()[i]->get_id()] = new_smallest->get_id(); //update previous node to node we removed
                            new_smallest->get_connections()[i]->set_tentative_distance(new_path_length); //update new shorter distance in min_pq
                            std::make_heap(const_cast<Node_VV**>(&min_pq.top()), const_cast<Node_VV**>(&min_pq.top())+ min_pq.size(), CompareDistance()); //sorts heap so that the closest by distance vertex is on top
                        }
                    }
                }
            }
        }
    }
    return previous_nodes;
}
void Path_search::multiple_points_path_search(std::vector<Node_VV *> user_stated_starting_points){
    result_json_file_alt << "\t\t\"paths\":[" << std::endl;
    int path_id = 1;
    std::ostringstream path_data;

    for(uint i = 0; i < user_stated_starting_points.size(); i++){
        std::map<int, int> tree_connections = simple_dijkstra_path_search(voronoi_diagram_vertices_, user_stated_starting_points[i]);
        for(uint j = 0; j < starting_points_.size(); j++){
            if(user_stated_starting_points[i]->get_id() != starting_points_[j]->get_id()){
                if(tree_connections.count(starting_points_[j]->get_id())){
                    path_data << "\t\t\t[";

                    int current = starting_points_[j]->get_id();
                    while(current != user_stated_starting_points[i]->get_id()){
                        int temp = tree_connections.find(current)->second;
                        Node_VV temp_node = indexed_nodes_.at(current);
                        path_data << temp_node.get_id() << ", ";
                        current = temp;
                    }
                    path_data << user_stated_starting_points[i]->get_id() << "]," << std::endl;
                    path_id++;
                }
            }
        }
    }

    std::string tmp_string = path_data.str();
    if(!tmp_string.empty()){
        tmp_string.pop_back();
    }
    if(!tmp_string.empty()){
        tmp_string.pop_back();
    }
    tmp_string.push_back('\n');

    result_json_file_alt << tmp_string;
    result_json_file_alt << "\t\t]" << std::endl;
}

void Path_search::multiple_points_path_search(){

    result_json_file_alt << "\t\t\"paths\":[" << std::endl;
    int path_id = 0;
    std::ostringstream path_data;

    for(uint i = 0; i < starting_points_.size(); i++){
        std::map<int, int> tree_connections = simple_dijkstra_path_search(voronoi_diagram_vertices_, starting_points_[i]);
        for(uint j = 0; j < starting_points_.size(); j++){
            if(starting_points_[i]->get_id() != starting_points_[j]->get_id()){
                if(tree_connections.count(starting_points_[j]->get_id())){
                    path_data << "\t\t\t[";

                    int current = starting_points_[j]->get_id();
                    while(current != starting_points_[i]->get_id()){
                        int temp = tree_connections.find(current)->second;
                        Node_VV temp_node = indexed_nodes_.at(current);
                        path_data << temp_node.get_id() << ", ";
                        current = temp;
                    }
                    path_data << starting_points_[i]->get_id() << "]," << std::endl;
                    path_id++;
                }
            }
        }
    }
    std::string tmp_string = path_data.str();
    tmp_string.pop_back();
    tmp_string.pop_back();
    tmp_string.push_back('\n');

    result_json_file_alt << tmp_string;
    result_json_file_alt << "\t\t]" << std::endl;
}

std::vector<Node_VV*> Path_search::get_starting_points(){
    return starting_points_;
}

std::vector<Node_VV*> Path_search::get_user_starting_points(){
    return user_starting_points_;
}

void Path_search::save_starting_points_to_pdb(std::string file_name){
    std::ofstream file;
    int counter = 0;
    file.open(file_name + ".pdb");
    for(int i = 0; i < starting_points_.size(); i++){
        file << "ATOM  " << std::setw(5) << starting_points_[i]->get_id() << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3) << starting_points_[i]->get_x() << std::setw(8) << std::setprecision(3) << starting_points_[i]->get_y() << std::setw(8) << std::setprecision(3) << starting_points_[i]->get_z() << "           H  " << std::endl;

    }
    file.close();
}

void Path_search::save_starting_points_to_pdb(std::string file_name, std::vector<Node_VV *> points_vector){
    std::ofstream file;
    int counter = 0;
    file.open(file_name + ".pdb");
    for(int i = 0; i < points_vector.size(); i++){
        file << "ATOM  " << std::setw(5) << points_vector[i]->get_id() << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3) << points_vector[i]->get_x() << std::setw(8) << std::setprecision(3) << points_vector[i]->get_y() << std::setw(8) << std::setprecision(3) << points_vector[i]->get_z() << "           H  " << std::endl;

    }
    file.close();
}


void Path_search::closest_common_starting_point(std::array<float, 3> qt_center){
    Node_VV* center_node_vv = NULL;
    float distance = std::numeric_limits<float>::infinity();
    for(int i = 0; i< voronoi_diagram_vertices_->size(); i++){
        if(!voronoi_diagram_vertices_->at(i).get_is_discarded()){
            float temp_distance = calculate_distance(voronoi_diagram_vertices_->at(i).get_coordinates(), qt_center);
            if(temp_distance < distance){
                distance = temp_distance;
                center_node_vv = &voronoi_diagram_vertices_->at(i);
            }
        }
    }

    simple_dijkstra_path_search(voronoi_diagram_vertices_, center_node_vv);

    for(int i = 0; i<voronoi_diagram_vertices_->size(); i++){
        if(!voronoi_diagram_vertices_->at(i).get_is_discarded()){
            if(voronoi_diagram_vertices_->at(i).get_is_starting_point()){
                for(int j = 0; j<voronoi_diagram_vertices_->at(i).get_connections().size(); j++){
                    if(!voronoi_diagram_vertices_->at(i).get_connections()[j]->get_is_discarded()){
                        if(voronoi_diagram_vertices_->at(i).get_connections()[j]->get_is_starting_point()){
                            if(voronoi_diagram_vertices_->at(i).get_connections()[j]->get_tentative_distance() < voronoi_diagram_vertices_->at(i).get_tentative_distance()){
                                voronoi_diagram_vertices_->at(i).set_is_starting_point(false);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Path_search::discard_starting_point_with_no_connections(){
    for(int i = 0; i<voronoi_diagram_vertices_->size(); i++){
        if(voronoi_diagram_vertices_->at(i).get_id() == 17934){
            Node_VV current = voronoi_diagram_vertices_->at(i);
        }
        bool has_valid_connections = false;
        if(voronoi_diagram_vertices_->at(i).get_is_starting_point()){
            for(int j = 0; j<voronoi_diagram_vertices_->at(i).get_connections().size(); j++){
                if(!voronoi_diagram_vertices_->at(i).get_connections()[j]->get_is_discarded()){
                    has_valid_connections = true;
                }
            }
            if(!has_valid_connections){
                voronoi_diagram_vertices_->at(i).set_is_starting_point(false);
            }
        }
    }
}

int Path_search::clusterize_starting_points(float epsilon, int minimal_no_points_in_cluster){
    int cluster_counter = 0;
    for(int i=0; i<starting_points_.size(); i++){
        if(starting_points_[i]->get_dbscan_label() != "undefined"){
            continue;
        }
        std::vector<Node_VV*> neighbours;
        for(uint j=0; j<starting_points_.size(); j++){
            if(starting_points_[i]->get_id() != starting_points_[j]->get_id()){
                if(calculate_distance(starting_points_[i]->get_coordinates(), starting_points_[j]->get_coordinates()) <= epsilon){
                    neighbours.push_back(starting_points_[j]);
                }
            }
        }

        if(neighbours.size() <= minimal_no_points_in_cluster){ //label as noise if points[i] doesn't have enough neighbours
            starting_points_[i]->set_dbscan_label("noise");
            continue;
        }
        cluster_counter++;
        starting_points_[i]->set_dbscan_label(std::to_string(cluster_counter));

        for(int j=0; j< neighbours.size(); j++){
            if(neighbours[j]->get_dbscan_label() == "noise"){
                neighbours[j]->set_dbscan_label(std::to_string(cluster_counter));
            }
            if(neighbours[j]->get_dbscan_label() != "undefined"){
                continue;
            }
            neighbours[j]->set_dbscan_label(std::to_string(cluster_counter));
            std::vector<Node_VV*> neighbours_to_add;
            for(uint k=0;k<starting_points_.size();k++){
                if(starting_points_[k]->get_id() != neighbours[j]->get_id()){
                    if(calculate_distance(starting_points_[k]->get_coordinates(), neighbours[j]->get_coordinates()) <= epsilon){
                        neighbours_to_add.push_back(starting_points_[k]);
                    }
                }
            }
            if(neighbours_to_add.size() >= minimal_no_points_in_cluster){
                neighbours.insert(neighbours.end(), neighbours_to_add.begin(), neighbours_to_add.end());
            }
        }
        for(int j=0; j<neighbours.size(); j++){
            neighbours[j]->set_is_starting_point(false);
        }
    }
    return cluster_counter;
}
