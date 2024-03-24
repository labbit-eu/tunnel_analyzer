#ifndef Point_3d_H
#define Point_3d_H

// class representing a single point for Voronoi diagram of spheres and Quasi-triangulation
#include <array>

class Point_3d{

    private:
        float x_, y_, z_;
        float radius_;

    public:
        Point_3d();
        Point_3d(float x, float y, float z);
        Point_3d(float x, float y, float z, float radius);

        float get_x();
        float get_y();
        float get_z();
        std::array<float,3> get_coordinates();
        float get_radius();

        void set_x(float x);
        void set_y(float y);
        void set_z(float z);
        void set_radius(float radius);
};

#endif