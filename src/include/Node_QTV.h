#ifndef Node_QTV_H
#define Node_QTV_H

#include "Point_3d.h"

class Node_QTV : public Point_3d{
    private:
        int qtv_id_;
        int qtv_atom_id_;

    public:
        Node_QTV();
        Node_QTV(int id, int atom_id, float x, float y, float z, float radius);

        int get_qtv_id();
        int get_qtv_atom_id();
        std::string print_to_pdb_string();
};

#endif