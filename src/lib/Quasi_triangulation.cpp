#include "../include/Quasi_triangulation.h"
#include <iterator>
#include <CGAL/iterator.h>
#include <array>
#include <vector>

void Quasi_triangulation::add_cgal_point(float x, float y, float z){
    cgal_quasi_triangulation_points_.push_back(Point(x,y,z));
}

void Quasi_triangulation::add_vertex(Node_QTV vertex){
    quasi_triangulation_vertices_.push_back(vertex);
}

void Quasi_triangulation::add_cell(Node_QTC cell){
    quasi_triangulation_cells_.push_back(cell);
}

void Quasi_triangulation::initialize_center_coordinates(std::array<float, 3> sums){
    center_coordinates_of_qt_ = {(sums[0]/quasi_triangulation_vertices_.size()), (sums[1]/quasi_triangulation_vertices_.size()), (sums[2]/quasi_triangulation_vertices_.size())};
}

void Quasi_triangulation::set_extreme_coordinates(float x, float neg_x, float y, float neg_y, float z, float neg_z){
    extreme_coordinates_ = {std::make_pair(x, neg_x), std::make_pair(y, neg_y), std::make_pair(z, neg_z)};
}

void Quasi_triangulation::save_alpha_shape_to_pdb(std::string file_name, Alpha_shape_3* alpha_shape){
    std::ofstream myfile_alpha;
    myfile_alpha.open (file_name + ".pdb");
    std::stringstream pts;
    std::stringstream ind;

    std::vector<Alpha_shape_3::Facet> facets;
    alpha_shape->get_alpha_shape_facets(std::back_inserter(facets), Alpha_shape_3::REGULAR);
    std::size_t nbf=facets.size();

    for (std::size_t i=0;i<nbf;++i){
        //To have a consistent orientation of the facet, always consider an exterior cell
        if ( alpha_shape->classify( facets[i].first )!=Alpha_shape_3::INTERIOR ){
            facets[i]=alpha_shape->mirror_facet( facets[i] );
        }

        int indices[3]={
            (facets[i].second+1)%4,
            (facets[i].second+2)%4,
            (facets[i].second+3)%4,
        };

            /// according to the encoding of vertex indices, this is needed to get
            /// a consistent orienation
        if ( facets[i].second%2==0 ) std::swap(indices[0], indices[1]);


        pts <<  "ATOM  " << std::setw(5) << (i*3) + 1 << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[0])->point().x() << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[0])->point().y() << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[0])->point().z() << "           H  " << "\n" <<
                "ATOM  " << std::setw(5) << (i*3) + 2 << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[1])->point().x() << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[1])->point().y() << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[1])->point().z() << "           H  " << "\n" <<
                "ATOM  " << std::setw(5) << (i*3) + 3 << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[2])->point().x() << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[2])->point().y() << std::setw(8) << std::setprecision(3) << facets[i].first->vertex(indices[2])->point().z() << "           H  " << "\n";

        ind <<  "CONECT" << std::setw(5) << 3*i+1 << std::setw(5) << 3*i+2 << "\n" <<
                "CONECT" << std::setw(5) << 3*i+1 << std::setw(5) << 3*i+3 << "\n" <<
                "CONECT" << std::setw(5) << 3*i+2 << std::setw(5) << 3*i+3 << "\n";
    }


    myfile_alpha << pts.str();
    myfile_alpha << ind.str();
    myfile_alpha.close();
}


void traverse(int current, std::vector<bool> &visited, int array_size, std::vector<std::vector <int>> &adj_matrix){
    visited[current] = true;
    for(int i = 0; i<array_size; i++){
        if(adj_matrix[current][i]){
            if(!visited[i]){
                traverse(i, visited, array_size, adj_matrix);
            }
        }
    }
}

bool has_hollow_alpha_shape(std::vector<Alpha_shape_facet> facets){
    int size = facets.size();
    std::vector<std::vector <int>> adjacency_matrix(size, std::vector<int>(size));
    std::vector<bool> visited(size);

    for(int i = 0; i< facets.size(); i++){
        for(int j = 0; j < facets.size(); j++){
            if(facets[i].check_if_neighbours(facets[j])){
                adjacency_matrix[i][j] = 1;
            }
        }
        visited[i] = false;
    }
    traverse(0, visited, facets.size(), adjacency_matrix);
    for(int i = 0; i<facets.size(); i++){
        if(visited[i] == false){
            return true;
        }
    }
    return false;
}

void Quasi_triangulation::calculate_alpha_shape(){
    //looks for an alpha_shape without hollows(scooped out insides) and creates concave hull facets based on cgal object information

    bool hull_has_hollows = false;
    int max_alpha = calculate_extreme_radius()*2;

    Alpha_shape_3 alpha_shape(cgal_quasi_triangulation_points_.begin(), cgal_quasi_triangulation_points_.end());

    while(!hull_has_hollows){
        alpha_shape.set_alpha(max_alpha);
        std::vector<Alpha_shape_3::Facet> cgal_facets;
        alpha_shape.get_alpha_shape_facets(std::back_inserter(cgal_facets), Alpha_shape_3::REGULAR);

        for(int i = 0; i < cgal_facets.size(); i++){
            if ( alpha_shape.classify(cgal_facets[i].first )!=Alpha_shape_3::EXTERIOR ){
                cgal_facets[i]=alpha_shape.mirror_facet(cgal_facets[i]);
            }
            CGAL_assertion(alpha_shape.classify( cgal_facets[i].first )==Alpha_shape_3::EXTERIOR  );
            int indices[3]={
                (cgal_facets[i].second+1)%4,
                (cgal_facets[i].second+2)%4,
                (cgal_facets[i].second+3)%4,
            };
            // according to the encoding of vertex indices, this is needed to get a consistent orienation
            if ( cgal_facets[i].second%2==0 ) std::swap(indices[0], indices[1]);
            std::array<float, 3> face_node_1{float(cgal_facets[i].first->vertex(indices[0])->point().x()), float(cgal_facets[i].first->vertex(indices[0])->point().y()), float(cgal_facets[i].first->vertex(indices[0])->point().z())};
            std::array<float, 3> face_node_2{float(cgal_facets[i].first->vertex(indices[1])->point().x()), float(cgal_facets[i].first->vertex(indices[1])->point().y()), float(cgal_facets[i].first->vertex(indices[1])->point().z())};
            std::array<float, 3> face_node_3{float(cgal_facets[i].first->vertex(indices[2])->point().x()), float(cgal_facets[i].first->vertex(indices[2])->point().y()), float(cgal_facets[i].first->vertex(indices[2])->point().z())};
            std::array<std::array<float, 3>, 3> my_facet{face_node_1, face_node_2, face_node_3};
            Alpha_shape_facet my_facet_object(i, my_facet);
            concave_hull_facets_.push_back(my_facet_object);
        }
        if(!has_hollow_alpha_shape(concave_hull_facets_)){
            concave_hull_facets_.clear();
            std::cout<<"alpha: "<<max_alpha<<std::endl;
            max_alpha = max_alpha - 3;
        }else{
            hull_has_hollows = true;
        }
    }
    while(hull_has_hollows){
        std::cout<<"alpha: "<<max_alpha<<std::endl;
        max_alpha++;
        alpha_shape.set_alpha(max_alpha);
        std::vector<Alpha_shape_3::Facet> cgal_facets;
        alpha_shape.get_alpha_shape_facets(std::back_inserter(cgal_facets), Alpha_shape_3::REGULAR);

        for(int i = 0; i < cgal_facets.size(); i++){
            if ( alpha_shape.classify(cgal_facets[i].first )!=Alpha_shape_3::EXTERIOR ){
                cgal_facets[i]=alpha_shape.mirror_facet(cgal_facets[i]);
            }
            CGAL_assertion(alpha_shape.classify( cgal_facets[i].first )==Alpha_shape_3::EXTERIOR  );
            int indices[3]={
                (cgal_facets[i].second+1)%4,
                (cgal_facets[i].second+2)%4,
                (cgal_facets[i].second+3)%4,
            };
            // according to the encoding of vertex indices, this is needed to get a consistent orienation
            if ( cgal_facets[i].second%2==0 ) std::swap(indices[0], indices[1]);
            std::array<float, 3> face_node_1{float(cgal_facets[i].first->vertex(indices[0])->point().x()), float(cgal_facets[i].first->vertex(indices[0])->point().y()), float(cgal_facets[i].first->vertex(indices[0])->point().z())};
            std::array<float, 3> face_node_2{float(cgal_facets[i].first->vertex(indices[1])->point().x()), float(cgal_facets[i].first->vertex(indices[1])->point().y()), float(cgal_facets[i].first->vertex(indices[1])->point().z())};
            std::array<float, 3> face_node_3{float(cgal_facets[i].first->vertex(indices[2])->point().x()), float(cgal_facets[i].first->vertex(indices[2])->point().y()), float(cgal_facets[i].first->vertex(indices[2])->point().z())};
            std::array<std::array<float, 3>, 3> my_facet{face_node_1, face_node_2, face_node_3};
            Alpha_shape_facet my_facet_object(i, my_facet);
            concave_hull_facets_.push_back(my_facet_object);
        }
        if(has_hollow_alpha_shape(concave_hull_facets_)){
            concave_hull_facets_.clear();
        }else{
            hull_has_hollows = false;
        }
    }
    alpha_value_ = max_alpha;
    alpha_shape.set_alpha(max_alpha + 10);
    std::vector<Alpha_shape_3::Facet> cgal_facets;
    alpha_shape.get_alpha_shape_facets(std::back_inserter(cgal_facets), Alpha_shape_3::REGULAR);
    for(int i = 0; i < cgal_facets.size(); i++){
        if ( alpha_shape.classify(cgal_facets[i].first )!=Alpha_shape_3::EXTERIOR ){
            cgal_facets[i]=alpha_shape.mirror_facet(cgal_facets[i]);
        }
        CGAL_assertion(alpha_shape.classify( cgal_facets[i].first )==Alpha_shape_3::EXTERIOR  );
        int indices[3]={
            (cgal_facets[i].second+1)%4,
            (cgal_facets[i].second+2)%4,
            (cgal_facets[i].second+3)%4,
        };
        // according to the encoding of vertex indices, this is needed to get a consistent orienation
        if ( cgal_facets[i].second%2==0 ) std::swap(indices[0], indices[1]);
        std::array<float, 3> face_node_1{float(cgal_facets[i].first->vertex(indices[0])->point().x()), float(cgal_facets[i].first->vertex(indices[0])->point().y()), float(cgal_facets[i].first->vertex(indices[0])->point().z())};
        std::array<float, 3> face_node_2{float(cgal_facets[i].first->vertex(indices[1])->point().x()), float(cgal_facets[i].first->vertex(indices[1])->point().y()), float(cgal_facets[i].first->vertex(indices[1])->point().z())};
        std::array<float, 3> face_node_3{float(cgal_facets[i].first->vertex(indices[2])->point().x()), float(cgal_facets[i].first->vertex(indices[2])->point().y()), float(cgal_facets[i].first->vertex(indices[2])->point().z())};
        std::array<std::array<float, 3>, 3> my_facet{face_node_1, face_node_2, face_node_3};
        Alpha_shape_facet my_facet_object(i, my_facet);
        larger_concave_hull_facets_.push_back(my_facet_object);
    }
}



std::vector<Node_QTV> Quasi_triangulation::get_qt_vertices(){
    return quasi_triangulation_vertices_;
}

std::vector<Node_QTC> Quasi_triangulation::get_qt_cells(){
    return quasi_triangulation_cells_;
}

std::array<float, 3> Quasi_triangulation::get_qt_center(){
    return center_coordinates_of_qt_;
}

std::vector<Alpha_shape_facet> Quasi_triangulation::get_concave_hull(){
    return concave_hull_facets_;
}

std::vector<Alpha_shape_facet> Quasi_triangulation::get_larger_concave_hull(){
    return larger_concave_hull_facets_;
}

std::array<std::pair<float, float>, 3> Quasi_triangulation::get_extreme_coordinates(){
    return extreme_coordinates_;
}

float Quasi_triangulation::calculate_extreme_radius(){
    float x_diameter = extreme_coordinates_[0].first - extreme_coordinates_[0].second;
    float y_diameter = extreme_coordinates_[1].first - extreme_coordinates_[1].second;
    float z_diameter = extreme_coordinates_[2].first - extreme_coordinates_[2].second;

    float current_max_diameter = std::max(x_diameter, y_diameter);

    return (std::max(current_max_diameter, z_diameter))/2;
}

void Quasi_triangulation::save_my_facets_to_pdb(){
    std::ofstream file;
    file.open("qt_facet_test.pdb");

    for(int i = 0; i< concave_hull_facets_.size();i++){
        file << concave_hull_facets_[i].print_to_pdb_string();
    }
    file.close();
}