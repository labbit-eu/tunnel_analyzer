#include "../include/Voronoi_diagram.h"
#include "../include/Utils.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
#include <queue>
#include <limits>
#include <algorithm>

#include <filesystem>
namespace fs = std::filesystem;
extern std::ofstream result_json_file_alt;



void Voronoi_diagram::add_vertex(Node_VV vertex){
    // adding single node to vertices vector
    voronoi_diagram_vertices_.push_back(vertex);
}

void Voronoi_diagram::save_diagram_to_pdb(std::string file_name, bool save_all){
    /*
    saving to pdb file for visualization purposes in pymol, if save_all = true (default) it saves every vertex,
    if save_all = false it will save only vertices that aren't discarded (Node_VV vertex.is_discarded = false)
    */
    std::ofstream file;
    std::stringstream nodes;
    std::stringstream connections;

    if(save_all){   // saving all vertices
        for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){    //  write node information in .pdb format
            nodes   << "ATOM  " << std::setw(5) << voronoi_diagram_vertices_[i].get_id()
                    << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                    << voronoi_diagram_vertices_[i].get_x() << std::setw(8) << std::setprecision(3)
                    << voronoi_diagram_vertices_[i].get_y() << std::setw(8) << std::setprecision(3)
                    << voronoi_diagram_vertices_[i].get_z() << "      " << std::setw(6) << std::setprecision(3)
                    << voronoi_diagram_vertices_[i].get_radius() << "           H  " << std::endl;
            connections << "CONECT" << std::setw(5) << voronoi_diagram_vertices_[i].get_id();
            for(uint j = 0; j < voronoi_diagram_vertices_[i].get_connections().size(); j++){    //  write node's connections information in .pdb format
                connections << std::setw(5) << voronoi_diagram_vertices_[i].get_connections()[j]->get_id();
            }
            connections << std::endl;
        }
    }else{  // saving only non discarded vertices
        for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){     //  write node information in .pdb format
            if(!voronoi_diagram_vertices_[i].get_is_discarded()){
                nodes   << "ATOM  " << std::setw(5) << voronoi_diagram_vertices_[i].get_id()
                        << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3)
                        << voronoi_diagram_vertices_[i].get_x() << std::setw(8) << std::setprecision(3)
                        << voronoi_diagram_vertices_[i].get_y() << std::setw(8) << std::setprecision(3)
                        << voronoi_diagram_vertices_[i].get_z() << "      " << std::setw(6) << std::setprecision(3)
                        << voronoi_diagram_vertices_[i].get_radius() << "           H  " << std::endl;
                connections << "CONECT" << std::setw(5) << voronoi_diagram_vertices_[i].get_id();
                for(uint j = 0; j < voronoi_diagram_vertices_[i].get_connections().size(); j++){    //  write node's connections information in .pdb format
                    if(!voronoi_diagram_vertices_[i].get_connections()[j]->get_is_discarded()){
                        connections << std::setw(5) << voronoi_diagram_vertices_[i].get_connections()[j]->get_id();
                    }
                }
                connections << std::endl;
            }
        }
    }

    file.open(file_name+".pdb");
    file << nodes.str();
    file << connections.str();
    file.close();
}

void Voronoi_diagram::save_diagram_to_mmcif_with_radius(std::string file_name, bool save_all){
    std::ofstream file;
    std::stringstream nodes;
    std::stringstream connections;

    nodes << "data_" << file_name << std::endl;
    nodes << "loop_\n";
    nodes << "_atom_site.id\n";
    nodes << "_atom_site.Cartn_x\n";
    nodes << "_atom_site.Cartn_y\n";
    nodes << "_atom_site.Cartn_z\n";
    nodes << "_atom_site.B_iso_or_equiv\n";

    if(save_all){   // saving all vertices
        for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){
            nodes << voronoi_diagram_vertices_[i].get_id() << " "
                << voronoi_diagram_vertices_[i].get_x() << " "
                << voronoi_diagram_vertices_[i].get_y() << " "
                << voronoi_diagram_vertices_[i].get_z() << " "
                << voronoi_diagram_vertices_[i].get_radius() << std::endl;
        }
    }else{  // saving only non discarded vertices
        for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){
            if(!voronoi_diagram_vertices_[i].get_is_discarded()){
                nodes << voronoi_diagram_vertices_[i].get_id() << " "
                      << voronoi_diagram_vertices_[i].get_x() << " "
                      << voronoi_diagram_vertices_[i].get_y() << " "
                      << voronoi_diagram_vertices_[i].get_z() << " "
                      << voronoi_diagram_vertices_[i].get_radius() << "\n";
            }
        }
    }
    file.open(file_name+".cif");
    file << nodes.str();
    file << connections.str();
    file.close();
}

void Voronoi_diagram::save_diagram_to_mmcif(std::string file_name, bool save_all){
    std::ofstream file;
    std::stringstream nodes;
    std::stringstream connections;

    nodes << "data_" << file_name << std::endl;
    nodes << "loop_\n";
    nodes << "_chem_comp_atom.atom_id\n";
    nodes << "_chem_comp_atom.model_Cartn_x\n";
    nodes << "_chem_comp_atom.model_Cartn_y\n";
    nodes << "_chem_comp_atom.model_Cartn_z\n";

    connections << "#\n";
    connections << "loop_\n";
    connections << "_chem_comp_bond.comp_id\n";
    connections << "_chem_comp_bond.atom_id_1\n";
    connections << "_chem_comp_bond.atom_id_2\n";

    if(save_all){   // saving all vertices
        for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){
            nodes << voronoi_diagram_vertices_[i].get_id() << " "
                << voronoi_diagram_vertices_[i].get_x() << " "
                << voronoi_diagram_vertices_[i].get_y() << " "
                << voronoi_diagram_vertices_[i].get_z() << std::endl;

            for(uint j = 0; j < voronoi_diagram_vertices_[i].get_connections().size(); j++){
                connections << "val " << voronoi_diagram_vertices_[i].get_id() << " " << voronoi_diagram_vertices_[i].get_connections()[j]->get_id() << std::endl;
            }
        }
    }else{  // saving only non discarded vertices
        for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){
            if(!voronoi_diagram_vertices_[i].get_is_discarded()){
                nodes << voronoi_diagram_vertices_[i].get_id() << " "
                      << voronoi_diagram_vertices_[i].get_x() << " "
                      << voronoi_diagram_vertices_[i].get_y() << " "
                      << voronoi_diagram_vertices_[i].get_z() << std::endl;
            }
            for(uint j = 0; j < voronoi_diagram_vertices_[i].get_connections().size(); j++){
                connections << "val " << voronoi_diagram_vertices_[i].get_id() << " " << voronoi_diagram_vertices_[i].get_connections()[j]->get_id() << std::endl;
            }
        }
    }
    file.open(file_name+".cif");
    file << nodes.str();
    file << connections.str();
    file.close();
}

void Voronoi_diagram::save_diagram_to_output_json(){
    std::ostringstream node_data;
    result_json_file_alt << "\t\t\"voronoi_nodes\":[" << std::endl;

    for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){
        if(!voronoi_diagram_vertices_[i].get_is_discarded()){
            node_data << "\t\t\t{" << std::endl;
                    node_data << "\t\t\t\t\"id\": " << voronoi_diagram_vertices_[i].get_id() << "," << std::endl;
                    node_data << "\t\t\t\t\"x\": " << voronoi_diagram_vertices_[i].get_x() << "," << std::endl;
                    node_data << "\t\t\t\t\"y\": " << voronoi_diagram_vertices_[i].get_y() << "," << std::endl;
                    node_data << "\t\t\t\t\"z\": " << voronoi_diagram_vertices_[i].get_z() << "," << std::endl;
                    node_data << "\t\t\t\t\"radius\": " << voronoi_diagram_vertices_[i].get_radius() << "," << std::endl;
                    node_data << "\t\t\t\t\"tangent_atoms\": " << "[";
                    for(uint j=0; j<voronoi_diagram_vertices_[i].get_tangent_qtvs().size(); j++){
                        node_data << voronoi_diagram_vertices_[i].get_tangent_qtvs()[j]->get_qtv_atom_id();
                        if(!((j+1) == voronoi_diagram_vertices_[i].get_tangent_qtvs().size())){
                            node_data << ", ";
                        }
                    }
                    node_data << "]" << std::endl;
            node_data << "\t\t\t}," << std::endl;
        }
    }
    std::string node_data_str = node_data.str();
    node_data_str.pop_back();
    node_data_str.pop_back();
    node_data_str.push_back('\n');

    result_json_file_alt << node_data_str;
    result_json_file_alt << "\t\t]," << std::endl;
}

void Voronoi_diagram::trim_discarded_vertices(){
    /*
    function removing discarded voronoi vertices from vertices vector
    and from other node's connections
    */
    for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){
        if(voronoi_diagram_vertices_[i].get_is_discarded()){
            for(uint j = 0; j < voronoi_diagram_vertices_[i].get_connections().size(); j++){
                voronoi_diagram_vertices_[i].get_connections()[j]->remove_connection(voronoi_diagram_vertices_[i].get_id());
            }
            voronoi_diagram_vertices_.erase(voronoi_diagram_vertices_.begin() + i);
            --i;
        }
    }
}

void Voronoi_diagram::filter_by_user_radius(float user_radius_filter_value){
    for(uint i = 0; i < voronoi_diagram_vertices_.size(); i++){
        if(voronoi_diagram_vertices_[i].get_radius() < user_radius_filter_value){
            voronoi_diagram_vertices_[i].set_is_discarded(true);
        }
    }
}

bool is_point_inside_triangle(float intersection[3], Alpha_shape_facet face){
    /*  check if given point lies within face boundaries

        input:
        intersection[3] -> point x,y,z of intersection with face's plane
        face -> face object defined by 3 points*/

    std::vector<float> u; //p2-p1
    std::vector<float> v; //p3-p1
    std::vector<float> n; //u x v -> cross product
    std::vector<float> w; //intersection - p1

    u = calculate_vector(face.get_vertices_of_face()[0][0], face.get_vertices_of_face()[0][1], face.get_vertices_of_face()[0][2],
                         face.get_vertices_of_face()[1][0], face.get_vertices_of_face()[1][1], face.get_vertices_of_face()[1][2]);
    v = calculate_vector(face.get_vertices_of_face()[0][0], face.get_vertices_of_face()[0][1], face.get_vertices_of_face()[0][2],
                         face.get_vertices_of_face()[2][0], face.get_vertices_of_face()[2][1], face.get_vertices_of_face()[2][2]);
    n = calculate_cross_product(u, v);
    w = calculate_vector(face.get_vertices_of_face()[0][0], face.get_vertices_of_face()[0][1], face.get_vertices_of_face()[0][2],
        intersection[0], intersection[1], intersection[2]);

    float gamma, beta, alfa;

    gamma = calculate_dot_product(calculate_cross_product(u, w), n) / calculate_dot_product(n, n);
    beta = calculate_dot_product(calculate_cross_product(w, v), n) / calculate_dot_product(n, n);
    alfa = 1 - gamma - beta;

    if((alfa<=1.0 && alfa>=0.0) && (beta<=1.0 && beta>=0.0) && (gamma<=1.0 && gamma>=0.0)){         //point inside
        return true;
    }else{                                                                                      //point outside
        return false;
    }
}

void Voronoi_diagram::filter_vertices_outside_of_concave_hull(std::vector<Alpha_shape_facet> concave_hull, std::array<float, 3> center){
    /*  Compare each vvertex with each convex hull face,
        check if intersection between segment (between vvertex and center of the protein) and plane (defined by hull's face) occurs,
        if it does, check if the intersection point lies within face surface boundaries, if it does increase counter++,
        after iterating through all faces check if counter%2 != 0, if true remove vertex.

        input:
        voronoi_diagram_vertices_ -> vector containing Vvertices
        concave_hull -> vector coontaining hull_faces
        center -> qtriangulation average center

        source: http://geomalgorithms.com/a05-_intersect-1.html*/
    std::vector<float> u;    //center - voronoi
    std::vector<float> w;    //voronoi - v[0] of plane
    float intersection_point[3];    //point of plane/segment intesection
    float D;
    float N;
    for(uint i = 0; i<voronoi_diagram_vertices_.size(); i++){
        if(!voronoi_diagram_vertices_[i].get_is_discarded()){                                                                       // for every vvertex
            u = calculate_vector(voronoi_diagram_vertices_[i].get_x(), voronoi_diagram_vertices_[i].get_y(), voronoi_diagram_vertices_[i].get_z(), center[0], center[1], center[2]);
            int counter = 0;
            for(uint j = 0; j<concave_hull.size(); j++){                                                                          // for every hull face
                w = calculate_vector(concave_hull[j].get_vertices_of_face()[0][0], concave_hull[j].get_vertices_of_face()[0][1], concave_hull[j].get_vertices_of_face()[0][2], voronoi_diagram_vertices_[i].get_x(), voronoi_diagram_vertices_[i].get_y(), voronoi_diagram_vertices_[i].get_z());
                D = calculate_dot_product(concave_hull[j].get_normal(), u);
                N = -1 *calculate_dot_product(concave_hull[j].get_normal(), w);

                if(fabs(D) < 0.000000001){        //segment is parallel to plane
                    if(N == 0){                       //segment lies in a plane
                        continue;
                    }else{
                        continue;                   //there is no intersection
                    }
                }
                float sI = N/D;
                if(sI<0 || sI>1){
                    continue;                       // no intersection
                }

                intersection_point[0] = sI * u[0] + voronoi_diagram_vertices_[i].get_x();                   // get intersection point coordinates
                intersection_point[1] = sI * u[1] + voronoi_diagram_vertices_[i].get_y();
                intersection_point[2] = sI * u[2] + voronoi_diagram_vertices_[i].get_z();

                if(is_point_inside_triangle(intersection_point,concave_hull[j])){                         //  check if point lies on face's surface
                    counter++;
                }
            }
            if(counter%2 != 0){                                                                         //  even counter means that the point lies outside convex hull -> remove this point 
                voronoi_diagram_vertices_[i].set_is_discarded(true);
                for(uint v=0; v<voronoi_diagram_vertices_[i].get_connections().size(); v++){
                    voronoi_diagram_vertices_[i].get_connections()[v]->set_is_starting_point(true);
                }
            }
        }
    }
}

void Voronoi_diagram::set_initial_voronoi_size(){
    voronoi_diagram_initial_size_ = voronoi_diagram_vertices_.size();
}

void Voronoi_diagram::calculate_starting_points_from_alpha_shapes(std::vector<Alpha_shape_facet> concave_hull, std::vector<Alpha_shape_facet> larger_concave_hull){
    std::set<std::array<std::array<float,3>,2>> edges;
    for(int i = 0; i< larger_concave_hull.size(); i++){
        for(int j = 0; j<3 ; j++){
            edges.insert(larger_concave_hull[i].get_edges()[j]);
        }
    }
    int initial_voronoi_size = voronoi_diagram_vertices_.size();
    int temp_id = initial_voronoi_size + 1;

    for(int i = 0; i< concave_hull.size(); i++){
        for(int j = 0; j<3 ; j++){
            if(edges.insert(concave_hull[i].get_edges()[j]).second){
                float avg_x = (concave_hull[i].get_edges()[j][0][0] + concave_hull[i].get_edges()[j][1][0]) / 2;
                float avg_y = (concave_hull[i].get_edges()[j][0][1] + concave_hull[i].get_edges()[j][1][1]) / 2;
                float avg_z = (concave_hull[i].get_edges()[j][0][2] + concave_hull[i].get_edges()[j][1][2]) / 2;
                Node_VV temp_node(temp_id, avg_x, avg_y, avg_z, 1);
                temp_id++;

                float closest_distance = std::numeric_limits<float>::infinity();;
                Node_VV *ptr;
                for(int k = 0; k <initial_voronoi_size; k++){
                    if(!voronoi_diagram_vertices_[k].get_is_discarded()){
                        float current_distance = calculate_distance(voronoi_diagram_vertices_[k].get_coordinates(), temp_node.get_coordinates());
                        if(current_distance < closest_distance){
                            closest_distance = current_distance;
                            ptr = &voronoi_diagram_vertices_[k];
                        }
                    }
                }
                temp_node.add_connection(ptr);
                temp_node.set_is_starting_point(true);
                temp_node.set_is_alpha_starting_point(true);
                ptr->add_connection(&temp_node);
                voronoi_diagram_vertices_.push_back(temp_node);
            }
        }
    }
}


void Voronoi_diagram::mark_user_starting_points(float search_radius, std::vector<std::array<float,3>> user_points){
    //this function sets to true additional bool flag on Node_VV.
    //I needed two seprate flags, one for user point and one for automatically detected starting points.
    for(int i=0; i<voronoi_diagram_vertices_.size(); i++){
        for(int j =0; j<user_points.size();j++){
            if(calculate_distance(voronoi_diagram_vertices_[i].get_coordinates(), user_points[j]) <= search_radius){
                voronoi_diagram_vertices_[i].set_is_user_starting_point(true);
            }
        }
    }
}

void Voronoi_diagram::search_for_user_starting_points(float search_radius, std::vector<std::array<float,3>> user_points){
    for(int i=0; i<voronoi_diagram_vertices_.size(); i++){
        voronoi_diagram_vertices_[i].set_is_starting_point(false);
        for(int j =0; j<user_points.size();j++){
            if(calculate_distance(voronoi_diagram_vertices_[i].get_coordinates(), user_points[j]) <= search_radius){
                voronoi_diagram_vertices_[i].set_is_starting_point(true);
            }
        }
    }
}

void normalize(std::vector<float> &vec){
    float mod = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    vec[0] = vec[0] / mod;
    vec[1] = vec[1] / mod;
    vec[2] = vec[2] / mod;
}

void calculate_rot_matrix(float rot_angle,std::vector<float> rot_vec, float rot_matrix[3][3]){
    float rot_cos = cos(rot_angle);
    float rot_sin = sin(rot_angle);
    float x = rot_vec[0];
    float y = rot_vec[1];
    float z = rot_vec[2];
    rot_matrix[0][0] = rot_cos + x*x*(1 - rot_cos);
    rot_matrix[0][1] = x*y*(1 - rot_cos) - z*rot_sin;
    rot_matrix[0][2] = y*rot_sin + x*z*(1 - rot_cos);
    rot_matrix[1][0] = z*rot_sin + x*y*(1 - rot_cos);
    rot_matrix[1][1] = rot_cos + y*y*(1 - rot_cos);
    rot_matrix[1][2] = -x*rot_sin + y*z*(1 - rot_cos);
    rot_matrix[2][0] = -y*rot_sin + x*z*(1 - rot_cos);
    rot_matrix[2][1] = x*rot_sin + y*z*(1 - rot_cos);
    rot_matrix[2][2] = rot_cos + z*z*(1 - rot_cos);
}

void rotate_point(float rot_mat[3][3], Node_QTV* in, Node_QTV &out){
    /*  rotates given point using previously calculated rotation matrix

        input:
        rot_mat -> rotation matrix
        in -> coordinates of a point to be rotated

        output:
        out -> coordinates of a rotated point
    */
    float x, y, z;
    x = rot_mat[0][0]*in->get_x() + rot_mat[0][1]*in->get_y() + rot_mat[0][2]*in->get_z();
    y = rot_mat[1][0]*in->get_x() + rot_mat[1][1]*in->get_y() + rot_mat[1][2]*in->get_z();
    z = rot_mat[2][0]*in->get_x() + rot_mat[2][1]*in->get_y() + rot_mat[2][2]*in->get_z();
    out = Node_QTV(in->get_qtv_id(), in->get_qtv_atom_id(), x, y, z, in->get_radius());
}

float operation_ACD(float var_x1, float var_x2, float var_y1, float var_y2, float var_r1, float var_r2){
    // operation used in calculate_delta()
    return pow(var_x1, 2) - pow(var_x2, 2) + pow(var_y1, 2) - pow(var_y2, 2) - pow(var_r1, 2) + pow(var_r2, 2); 
}

void calculate_delta(Node_QTV v1, Node_QTV v2, Node_QTV v3, float &x, float &y,float &z,float &r,bool &is_ok){  //  Apollonius problem
    /*  Function solving set of equations of Apollonius Problem(find circle tangent to three other circles):
        (xs - x1)^2 + (ys - y1)^2 = (rs - s1r1)^2
        (xs - x2)^2 + (ys - y2)^2 = (rs - s2r2)^2
        (xs - x3)^2 + (ys - y3)^2 = (rs - s3r3)^2

        xs, ys, rs -> solution values

        input:
        v1, v2, v3 -> QTnodes of shared QTface(gate)

        output:
        x, y, z, r -> values of a minima
        is_ok -> boolean to check if solving equations was possible(quadratic equation)
    */
    float a = v1.get_x() - v2.get_x();
    float b = v1.get_x() - v3.get_x();
    float c = v1.get_y() - v2.get_y();
    float d = v1.get_y() - v3.get_y();

    float e = (-1*v1.get_radius()) - (-1*v2.get_radius());
    float f = (-1*v1.get_radius()) - (-1*v3.get_radius());

    float m = c - ((d*a) / b);

    //  creating auxiliary variables
    float A = operation_ACD(v1.get_x(), v2.get_x(), v1.get_y(), v2.get_y(), v1.get_radius(), v2.get_radius());
    float B = A * -1;
    float C = operation_ACD(v1.get_x(), v3.get_x(), v1.get_y(), v3.get_y(), v1.get_radius(), v3.get_radius());
    float E = C * -1;

    float Q = (e / m) - ((f*a) / (b*m));
    float N = (A / (2*m)) + ((a*E) / (2*b*m));
    float P = (e / a) - ((c*e) / (a*m)) + ((c*f) / (b*m));
    float M = (A / (2*a)) + ((c*B) / (2*a*m)) + ((c*C) / (2*b*m));

    float delta_a = pow(P, 2) + pow(Q, 2) - 1;
    float delta_b = (2*M*P) - (2*v1.get_x()*P) + (2*N*Q) - (2*v1.get_y()*Q) + (2*(-1)*v1.get_radius());
    float delta_c = pow(M, 2) - (2*v1.get_x()*M) + pow(v1.get_x(), 2) + pow(N, 2) - (2*v1.get_y()*N) + pow(v1.get_y(), 2) - pow(v1.get_radius(), 2);

    float delta = pow(delta_b, 2) - (4*delta_a*delta_c);

    float result1 = ((-delta_b - sqrt(delta)) / (2*delta_a));       //  two possible results
    float result2 = ((-delta_b + sqrt(delta)) / (2*delta_a));

    if(result1>0 && result2>0){         //always take smaller result -> circle tangent externally -> smallest radius possible
        if(result1 < result2){
            x = M + P*result1;
            y = N + Q*result1;
            z = v1.get_z();
            r = result1;
            is_ok = true;
        }else{
            x = M + P*result2;
            y = N + Q*result2;
            z = v1.get_z();
            r = result2;
            is_ok = true;
        }
    }else if(result1>0 && result2<0){
        x = M + P*result1;
        y = N + Q*result1;
        z = v1.get_z();
        r = result1;
        is_ok = true;
    }else if(result2>0 && result1<0){
        x = M + P*result2;
        y = N + Q*result2;
        z = v1.get_z();
        r = result2;
        is_ok = true;
    }else{                      // no solution if result1 < 0 && result2 < 0
        is_ok = false;
    }

}

void rotate_back(float &x, float &y, float &z,float rot_angle, std::vector<float> rot_vec){
    /*  Function to rotate back newly calculated minima's values

        input:
        x,y,z -> coordinates of a point to be rotated
        rot_angle -> rotation angle for calculating reversed rotation matrix
        rot_vec -> rotation vector for calculating reversed rotation matrix

        output:
        x,y,z -> valus of a point that was rotated back
    */
    float rot_matrix[3][3];
    float x_temp, y_temp, z_temp;
    x_temp = x;
    y_temp = y;
    z_temp = z;
    rot_vec[0] = -rot_vec[0];
    rot_vec[1] = -rot_vec[1];
    rot_vec[2] = -rot_vec[2];
    calculate_rot_matrix(rot_angle, rot_vec, rot_matrix);

    x = rot_matrix[0][0]*x_temp + rot_matrix[0][1]*y_temp + rot_matrix[0][2]*z_temp;
    y = rot_matrix[1][0]*x_temp + rot_matrix[1][1]*y_temp + rot_matrix[1][2]*z_temp;
    z = rot_matrix[2][0]*x_temp + rot_matrix[2][1]*y_temp + rot_matrix[2][2]*z_temp;
}

void rotation(Node_QTV* s1, Node_QTV* s2, Node_QTV* s3, float &x, float &y, float &z, float &r, bool &is_ok){
    std::vector<float> ab, ac;
    ab = calculate_vector(s1->get_x(), s1->get_y(), s1->get_z(), s2->get_x(), s2->get_y(), s2->get_z());
    ac = calculate_vector(s1->get_x(), s1->get_y(), s1->get_z(), s3->get_x(), s3->get_y(), s3->get_z());

    std::vector<float> face_normal = calculate_cross_product(ab, ac);
    normalize(face_normal);
    std::vector<float> rot_vec = calculate_cross_product(face_normal, {0, 0, 1});
    normalize(rot_vec);

    float rot_angle = acos(calculate_dot_product(face_normal, {0, 0, 1}));
    float rot_matrix[3][3];
    calculate_rot_matrix(rot_angle, rot_vec, rot_matrix);   //  rotation matrix needed to rotate the face in such a way that z-coordinate is the same for each point

    Node_QTV rotated_s1, rotated_s2, rotated_s3;
    rotate_point(rot_matrix, s1, rotated_s1);
    rotate_point(rot_matrix, s2, rotated_s2);
    rotate_point(rot_matrix, s3, rotated_s3);

    calculate_delta(rotated_s1, rotated_s2, rotated_s3, x, y, z, r, is_ok);
    if(is_ok){  //  only rotate back if calculate_delta() ended with a correct result
        rotate_back(x, y, z, rot_angle, rot_vec);
    }

}

void Voronoi_diagram::calculate_minimas(int user_radius_filter_value){
    int voronoi_size = voronoi_diagram_vertices_.size();
    int temp_id = voronoi_size + 1;
    for(int i=0; i<voronoi_size;i++){
        if(!voronoi_diagram_vertices_[i].get_is_discarded()){
            for(int j = 0; j< voronoi_diagram_vertices_[i].get_connections().size();j++){
                if(!voronoi_diagram_vertices_[i].get_connections()[j]->get_is_discarded()){
                    if(!voronoi_diagram_vertices_[i].get_connections()[j]->get_is_minima()){
                        float x = 0, y = 0, z = 0, r = 0;
                        bool is_ok;
                        rotation(voronoi_diagram_vertices_[i].get_gates()[j].get_gate_qtv_nodes()[0],
                                voronoi_diagram_vertices_[i].get_gates()[j].get_gate_qtv_nodes()[1],
                                voronoi_diagram_vertices_[i].get_gates()[j].get_gate_qtv_nodes()[2],
                                x, y, z, r, is_ok);

                        if((x > voronoi_diagram_vertices_[i].get_x() && x > voronoi_diagram_vertices_[i].get_connections()[j]->get_x()) || (x < voronoi_diagram_vertices_[i].get_x() && x < voronoi_diagram_vertices_[i].get_connections()[j]->get_x())){
                            is_ok = false;
                        }

                        if(is_ok){
                            Node_VV temp_node(temp_id, x, y, z, r, true);
                            temp_node.add_connection(&voronoi_diagram_vertices_[i]);
                            temp_node.add_connection(voronoi_diagram_vertices_[i].get_connections()[j]);
                            temp_node.set_tangent_qtv_at(voronoi_diagram_vertices_[i].get_gates()[j].get_gate_qtv_nodes()[0]);
                            temp_node.set_tangent_qtv_at(voronoi_diagram_vertices_[i].get_gates()[j].get_gate_qtv_nodes()[1]);
                            temp_node.set_tangent_qtv_at(voronoi_diagram_vertices_[i].get_gates()[j].get_gate_qtv_nodes()[2]);

                            if(temp_node.get_radius() < user_radius_filter_value){
                                temp_node.set_is_discarded(true);
                            }
                            if(voronoi_diagram_vertices_[i].get_is_starting_point() && voronoi_diagram_vertices_[i].get_connections()[j]->get_is_starting_point()){
                                temp_node.set_is_starting_point(true);
                            }

                            voronoi_diagram_vertices_.push_back(temp_node);

                            voronoi_diagram_vertices_[i].get_connections()[j]->swap_connection(voronoi_diagram_vertices_[i].get_id(),&voronoi_diagram_vertices_.back());
                            voronoi_diagram_vertices_[i].swap_connection(voronoi_diagram_vertices_[i].get_connections()[j]->get_id(), &voronoi_diagram_vertices_.back());

                            temp_id++;
                        }
                    }
                }
            }
        }
    }
}

void Voronoi_diagram::calculate_distances_to_neighbours(){
    for(int i = 0; i< voronoi_diagram_vertices_.size(); i++){
        voronoi_diagram_vertices_[i].calculate_distance_to_neighbours();
    }
}

void Voronoi_diagram::save_starting_points_to_pdb(std::string file_name){
    std::ofstream file;
    int counter = 0;
    file.open(file_name + ".pdb");
    for(int i = 0; i < voronoi_diagram_vertices_.size(); i++){
        if(!voronoi_diagram_vertices_[i].get_is_discarded()){
            if(voronoi_diagram_vertices_[i].get_is_starting_point()){
                file << "ATOM  " << std::setw(5) << voronoi_diagram_vertices_[i].get_id() << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3) << voronoi_diagram_vertices_[i].get_x() << std::setw(8) << std::setprecision(3) << voronoi_diagram_vertices_[i].get_y() << std::setw(8) << std::setprecision(3) << voronoi_diagram_vertices_[i].get_z() << "           H  " << std::endl;
                counter ++;
            }
        }
    }
    file.close();
}

void Voronoi_diagram::save_minimas_to_pdb(std::string file_name){
    std::ofstream file;
    int counter = 0;
    file.open(file_name + ".pdb");
    for(int i = 0; i < voronoi_diagram_vertices_.size(); i++){
        if(voronoi_diagram_vertices_[i].get_is_minima()){
            file << "ATOM  " << std::setw(5) << voronoi_diagram_vertices_[i].get_id() << "  H   UNK A   1    " << std::fixed << std::setw(8) << std::setprecision(3) << voronoi_diagram_vertices_[i].get_x() << std::setw(8) << std::setprecision(3) << voronoi_diagram_vertices_[i].get_y() << std::setw(8) << std::setprecision(3) << voronoi_diagram_vertices_[i].get_z() << "           H  " << std::endl;
            counter ++;
        }
    }
    file.close();
}

std::vector<Node_VV> Voronoi_diagram::get_voronoi_diagram(){
    return voronoi_diagram_vertices_;
}

std::vector<Node_VV> *Voronoi_diagram::get_voronoi_diagram_adress(){
    return &voronoi_diagram_vertices_;
}
