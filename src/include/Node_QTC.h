#ifndef NODE_QTC
#define NODE_QTC

class Node_QTC{
    private:
        int qtc_id_;
        int vv_ids_[4];
        int qtc_ids_[4];

    public:
        Node_QTC();
        Node_QTC(int id);
        Node_QTC(int id, int vv_id1, int vv_id2, int vv_id3, int vv_id4, int qtc_id1, int qtc_id2, int qtc_id3, int qtc_id4);

        int id();
        int* get_vv_ids();
        int* get_qtc_ids();
};

#endif