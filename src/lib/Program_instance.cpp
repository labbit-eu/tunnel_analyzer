#include "../include/Program_instance.h"

#include "../include/Utils.h"
#include <fstream>
#include <limits>

Instance::Instance(std::string snapshot_qtf_file_name){
    qtf_file_name_ = snapshot_qtf_file_name;
    parse_user_input_file("input.txt");

    parse_qtf_file();
    set_connections();

    std::cout<<"USER RADIUS: " << radius_filter_value_ << std::endl;
    voronoi_diagram_.filter_by_user_radius(radius_filter_value_);

    std::cout<<"Alpha parameter estimation..."<<std::endl;
    quasi_triangulation_.calculate_alpha_shape();

    std::cout<<"Out of hull filtering..."<<std::endl;
    voronoi_diagram_.filter_vertices_outside_of_concave_hull(quasi_triangulation_.get_concave_hull(), quasi_triangulation_.get_qt_center());

    voronoi_diagram_.calculate_minimas(radius_filter_value_);

    voronoi_diagram_.save_diagram_to_output_json();

    std::cout<<"Calculating distances..."<<std::endl;
    voronoi_diagram_.calculate_distances_to_neighbours();

    if(user_points_of_search.empty()){      //automatic search
        automatic_search();
    }else if(user_points_of_search.size() == 1){
        single_starting_sphere_search();
    }else{
        multiple_starting_spheres_search();
    }
}

void Instance::automatic_search(){
    std::cout << "Automatic starting point detection." << std::endl;

    Path_search paths(voronoi_diagram_.get_voronoi_diagram_adress());
    paths.discard_starting_point_with_no_connections();
    paths.closest_common_starting_point(quasi_triangulation_.get_qt_center());
    paths.save_starting_points_to_vector();
    std::cout<<"Size of starting points: "<< paths.get_starting_points().size() << std::endl;

    int no_clusters = paths.clusterize_starting_points(2*radius_filter_value_, 2);
    paths.save_starting_points_to_vector();

    std::cout<<"CLUSTERS: "<<no_clusters<<std::endl;
    std::cout<<"Size of starting points after clustering: "<< paths.get_starting_points().size() << std::endl;
    std::cout<<"Searching for paths within snapshot."<<std::endl;

    paths.multiple_points_path_search();
    std::cout<<"Search done."<<std::endl;
}

void Instance::single_starting_sphere_search(){
    std::cout << "One starting sphere stated - searching for possible starting points." << std::endl;

    Path_search paths(voronoi_diagram_.get_voronoi_diagram_adress());
    paths.discard_starting_point_with_no_connections();
    paths.closest_common_starting_point(quasi_triangulation_.get_qt_center());
    paths.save_starting_points_to_vector();

    int no_clusters = paths.clusterize_starting_points(2*radius_filter_value_, 2);
    paths.save_starting_points_to_vector();

    std::cout<<"Size of starting points after clustering: "<< paths.get_starting_points().size() << std::endl;
    std::cout<<"Size of starting points: "<< paths.get_starting_points().size() << std::endl;

    voronoi_diagram_.mark_user_starting_points(4, user_points_of_search);
    paths.save_starting_points_to_vector();

    std::cout<<"Size of user starting points: "<< paths.get_user_starting_points().size() << std::endl;
    std::cout<<"Searching for paths within snapshot."<<std::endl;
    paths.multiple_points_path_search(paths.get_user_starting_points());
    std::cout<<"Search done."<<std::endl;
}

void Instance::multiple_starting_spheres_search(){
    std::cout << "More than one starting sphere stated - searching for possible starting points." << std::endl;
    Path_search paths(voronoi_diagram_.get_voronoi_diagram_adress());
    voronoi_diagram_.search_for_user_starting_points(2, user_points_of_search);
    paths.save_starting_points_to_vector();
    std::cout<<"Size of starting points: "<< paths.get_starting_points().size() << std::endl;

    int no_clusters = paths.clusterize_starting_points(radius_filter_value_, 2);
    paths.save_starting_points_to_vector();

    std::cout<<"Size of starting points after clustering: "<< paths.get_starting_points().size() << std::endl;

    std::cout<<"Searching for paths within snapshot."<<std::endl;
    paths.multiple_points_path_search();
    std::cout<<"Search done."<<std::endl;
}

void Instance::parse_user_input_file(std::string user_input_file_name){
    /*
    Parsing function for user input file.
    Sets qtf file name and filter radius value.
    */

    std::ifstream input_file;
    input_file.open(user_input_file_name);

    if(input_file.is_open()){
        std::string single_file_line;
        getline(input_file, single_file_line);
        getline(input_file, single_file_line);
        radius_filter_value_ = stof(single_file_line);
        while(getline(input_file, single_file_line)){
            std::istringstream coords(single_file_line);
            float x,y,z;
            if(!(coords >> x >> y >> z)){
                break;
            }
            user_points_of_search.push_back({x, y, z});
            std::cout << "X: " << x << " Y: " << y << " Z: " << z << std::endl;
        }
    }else{
        std::cout << "User input file not found!" << std::endl;
    }

    input_file.close();
}

void Instance::parse_qtf_file(){
    // parsing function for .qtf file, also sets average quasi triangulation center
    float x_sum = 0, y_sum = 0, z_sum = 0;
    float extreme_x = std::numeric_limits<float>::min(),
        extreme_neg_x = std::numeric_limits<float>::max(),
        extreme_y = std::numeric_limits<float>::min(),
        extreme_neg_y = std::numeric_limits<float>::max(),
        extreme_z = std::numeric_limits<float>::min(),
        extreme_neg_z = std::numeric_limits<float>::max();
    std::cout << "QTF FILE NAME: " << qtf_file_name_ << std::endl;
    std::ifstream qtf_file(qtf_file_name_);
    if (qtf_file.is_open()){
        std::string temp_line, QTVTX = "QTVTX", QTCELL = "QTCELL", VDVTX = "VDVTX";
        std::vector<std::string> chunks;
        while (getline(qtf_file, temp_line)){
            chunks = split(temp_line, '\t');
            if (chunks[0].compare(QTVTX) == 0){     // save QTvertices
                uint qtv_id, qtv_atid;
                float qtv_x, qtv_y, qtv_z, qtv_rad;
                std::stringstream s_qtv_id(chunks[1]);
                std::stringstream s_qtv_atid(chunks[7]);
                std::stringstream s_qtv_x(chunks[3]);
                std::stringstream s_qtv_y(chunks[4]);
                std::stringstream s_qtv_z(chunks[5]);
                std::stringstream s_qtv_rad(chunks[6]);
                s_qtv_id >> qtv_id;
                s_qtv_atid >> qtv_atid;
                s_qtv_x >> qtv_x;
                s_qtv_y >> qtv_y;
                s_qtv_z >> qtv_z;
                s_qtv_rad >> qtv_rad;
                x_sum += qtv_x;
                y_sum += qtv_y;
                z_sum += qtv_z;
                Node_QTV temp_node(qtv_id, qtv_atid, qtv_x, qtv_y, qtv_z, qtv_rad);
                quasi_triangulation_.add_vertex(temp_node);
                quasi_triangulation_.add_cgal_point(qtv_x,qtv_y, qtv_z);
                if(qtv_x > extreme_x){
                    extreme_x = qtv_x;
                }
                if(qtv_x < extreme_neg_x){
                    extreme_neg_x = qtv_x;
                }
                if(qtv_y > extreme_y){
                    extreme_y = qtv_y;
                }
                if(qtv_y < extreme_neg_y){
                    extreme_neg_y = qtv_y;
                }
                if(qtv_z > extreme_z){
                    extreme_z = qtv_z;
                }
                if(qtv_z < extreme_neg_z){
                    extreme_neg_z = qtv_z;
                }

            } else if (chunks[0].compare(QTCELL) == 0){   //  save QTcells
                uint qtc_id, vv_id1, vv_id2, vv_id3, vv_id4;
                uint qtc_id1, qtc_id2, qtc_id3, qtc_id4;
                std::stringstream s_qtc_id(chunks[1]);
                std::stringstream s_vv_id1(chunks[2]);
                std::stringstream s_vv_id2(chunks[3]);
                std::stringstream s_vv_id3(chunks[4]);
                std::stringstream s_vv_id4(chunks[5]);
                std::stringstream s_qtc_id1(chunks[6]);
                std::stringstream s_qtc_id2(chunks[7]);
                std::stringstream s_qtc_id3(chunks[8]);
                std::stringstream s_qtc_id4(chunks[9]);
                s_qtc_id >> qtc_id;
                s_vv_id1 >> vv_id1;
                s_vv_id2 >> vv_id2;
                s_vv_id3 >> vv_id3;
                s_vv_id4 >> vv_id4;
                s_qtc_id1 >> qtc_id1;
                s_qtc_id2 >> qtc_id2;
                s_qtc_id3 >> qtc_id3;
                s_qtc_id4 >> qtc_id4;
                Node_QTC temp_node(qtc_id, vv_id1, vv_id2, vv_id3, vv_id4, qtc_id1, qtc_id2, qtc_id3, qtc_id4);
                quasi_triangulation_.add_cell(temp_node);
            } else if (chunks[0].compare(VDVTX) == 0){    //  save voronoi vertices
                uint vv_id;
                float vv_x, vv_y, vv_z;
                std::stringstream s_vv_id(chunks[1]);
                std::stringstream s_vv_x(chunks[2]);
                std::stringstream s_vv_y(chunks[3]);
                std::stringstream s_vv_z(chunks[4]);
                s_vv_id >> vv_id;
                s_vv_x >> vv_x;
                s_vv_y >> vv_y;
                s_vv_z >> vv_z;
                Node_VV temp_node(vv_id, vv_x, vv_y, vv_z);
                voronoi_diagram_.add_vertex(temp_node);
            }
        }
        quasi_triangulation_.initialize_center_coordinates({x_sum, y_sum, z_sum});
        quasi_triangulation_.set_extreme_coordinates(extreme_x, extreme_neg_x, extreme_y, extreme_neg_y, extreme_z, extreme_neg_z);
    }else{
        std::cout << "QTF file not found!" << std::endl;
    }
}

void Instance::set_connections(){
    /*
    Setting connections between voronoi diagram vertices,
    iterate through qt-cells to add the ids of the atoms tangent to each
    vv-node, also add the connections to the vv-nodes
    */
    for(uint i = 0; i < voronoi_diagram_.voronoi_diagram_vertices_.size();i++){
        for(short j = 0; j < 4; j++){
              voronoi_diagram_.voronoi_diagram_vertices_[i].set_tangent_qtv_at(&quasi_triangulation_.quasi_triangulation_vertices_[quasi_triangulation_.quasi_triangulation_cells_[i].get_vv_ids()[j] - 1]);
        }
        /*
        After adding the quasi-triangulation vertices, create the gates
        of distances from the vv node to each quasi-triangulation vertex
        minus the radius of the quasi-triangulation vertex
        */
        voronoi_diagram_.voronoi_diagram_vertices_[i].update_radius();
        voronoi_diagram_.voronoi_diagram_vertices_[i].create_gates();
        for(short j = 0; j < 4; j++){
            if(quasi_triangulation_.quasi_triangulation_cells_[i].get_qtc_ids()[j] > 0){
                voronoi_diagram_.voronoi_diagram_vertices_[i].add_connection(&voronoi_diagram_.voronoi_diagram_vertices_[quasi_triangulation_.quasi_triangulation_cells_[i].get_qtc_ids()[j] - 1]);
                voronoi_diagram_.voronoi_diagram_vertices_[i].get_gates()[j].set_end_vv_node(&voronoi_diagram_.voronoi_diagram_vertices_[quasi_triangulation_.quasi_triangulation_cells_[i].get_qtc_ids()[j] - 1]);
            }
        }
        /*
        Initial radius check, due to possible errors while generating .qtf file
        some of the vertices might have negative radius, mark vertex with negative
        radius (bool is_discarded = true)
        */
        if(voronoi_diagram_.voronoi_diagram_vertices_[i].get_radius() < 0 ){
            voronoi_diagram_.voronoi_diagram_vertices_[i].set_is_discarded(true);
        }
        /*
        Initial out of bounds filtering using extreme coordinates of quasi-triangulation
        (bool is_discarded = true)
        */
        std::array<float,3> current_coordinates = voronoi_diagram_.voronoi_diagram_vertices_[i].get_coordinates();
        if(current_coordinates[0] > quasi_triangulation_.extreme_coordinates_[0].first){
            voronoi_diagram_.voronoi_diagram_vertices_[i].set_is_discarded(true);
        }else if(current_coordinates[0] < quasi_triangulation_.extreme_coordinates_[0].second){
            voronoi_diagram_.voronoi_diagram_vertices_[i].set_is_discarded(true);
        }else if(current_coordinates[1] > quasi_triangulation_.extreme_coordinates_[1].first){
            voronoi_diagram_.voronoi_diagram_vertices_[i].set_is_discarded(true);
        }else if(current_coordinates[1] < quasi_triangulation_.extreme_coordinates_[1].second){
            voronoi_diagram_.voronoi_diagram_vertices_[i].set_is_discarded(true);
        }else if(current_coordinates[2] > quasi_triangulation_.extreme_coordinates_[2].first){
            voronoi_diagram_.voronoi_diagram_vertices_[i].set_is_discarded(true);
        }else if(current_coordinates[2] < quasi_triangulation_.extreme_coordinates_[2].second){
            voronoi_diagram_.voronoi_diagram_vertices_[i].set_is_discarded(true);
        }
    }
}
