#include "../include/Node_QTC.h"

Node_QTC::Node_QTC(){
    qtc_id_ = -1;
}

Node_QTC::Node_QTC(int id){
    qtc_id_ = id;
}

Node_QTC::Node_QTC(int id, int vv_id1, int vv_id2, int vv_id3, int vv_id4, int qtc_id1, int qtc_id2, int qtc_id3, int qtc_id4){
    qtc_id_ = id;
    vv_ids_[0] = vv_id1;
    vv_ids_[1] = vv_id2;
    vv_ids_[2] = vv_id3;
    vv_ids_[3] = vv_id4;
    qtc_ids_[0] = qtc_id1;
    qtc_ids_[1] = qtc_id2;
    qtc_ids_[2] = qtc_id3;
    qtc_ids_[3] = qtc_id4;
}

int Node_QTC::id(){
    return qtc_id_;
}

int* Node_QTC::get_vv_ids(){
    return vv_ids_;
}

int* Node_QTC::get_qtc_ids(){
    return qtc_ids_;
}