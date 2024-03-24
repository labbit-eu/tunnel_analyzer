#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>

#include "src/include/Utils.h"
#include "src/include/Program_instance.h"

// Automated method for analysis of dynamic pathways in proteins
std::ofstream result_json_file_alt;

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
    result_json_file_alt.open("output_alt.json");
    result_json_file_alt << "{" << std::endl;
    result_json_file_alt << "\t\"snapshots\": [" << std::endl;
    std::cout << "Tunnel Analyzer" << std::endl;

    std::ifstream input_file;
    std::string directory_path;
    input_file.open("input.txt");
    if(input_file.is_open()){
        std::string single_file_line;
        getline(input_file, single_file_line);
        directory_path = single_file_line;
    }

    if(!fs::exists(directory_path) or !fs::is_directory(directory_path)) {
        std::cout << "Path is not a directory or does not exist." << std::endl;
        return 0;
    }

    int file_number = 1;
    size_t count = 0;
    for (const auto &entry : fs::directory_iterator(directory_path)) {
        if (fs::is_regular_file(entry.path())) {
            result_json_file_alt << "\t\t{" << std::endl;
            result_json_file_alt << "\t\t\"snapshot_id\": " << file_number << "," << std::endl;

            std::cout << "Run: " << file_number << std::endl;
            std::string qtf_file = directory_path + "/" +entry.path().filename().string();
            Instance instance(qtf_file);
            count++;
            if(count != std::distance(fs::directory_iterator(directory_path), fs::directory_iterator{})){
                result_json_file_alt << "\t\t}," << std::endl;
            }else{
                result_json_file_alt << "\t\t}" << std::endl;
            }
            file_number ++;
        }
    }
    std::cout << "Finished" << std::endl;

    result_json_file_alt << "\t]" << std::endl;
    result_json_file_alt << "}" << std::endl;
    result_json_file_alt.close();
}
