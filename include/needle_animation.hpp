#ifndef NEEDLE_ANIMATION_H
#define NEEDLE_ANIMATION_H

#include <iostream>
#include <chrono>
#include <armadillo>
#include <vector>

#include "euler_rotations.hpp"
#include "gnuplot-iostream.h"
#include "needle_geometry.hpp"

class NeedleAnimation: public NeedleGeometry
{

public:
    NeedleAnimation();

    void animate(arma::fvec roa_f_f, arma::fvec euler_angles, 
        std::vector<std::tuple<float, float, float>> beam_points);

private:

    Gnuplot m_gp;
    
    void create_handle(arma::fvec center, arma::fvec orientation);
    void create_needle(std::vector<std::vector<arma::fvec>> beam_points);
    
    
    std::vector<arma::fmat> cylinder(arma::fvec center, arma::fvec orientation);
    std::vector<std::vector<std::string>> polygons_to_strings(std::vector<arma::fmat> polygons);
    void polygon_plot(std::vector<std::vector<std::string>> str_polygon);

};


NeedleAnimation::NeedleAnimation()
{
    m_gp << "set style line 1 lc rgb 'black' lw 1.5\n";
    m_gp << "set xrange [-0.4:0.4]\nset yrange [-0.4:0.4]\nset zrange [-0.4:0.4]\n";
    m_gp << "set hidden3d\n";
    m_gp << "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\n";
    m_gp << "set view 90, 0, 1, 1\n";
    m_gp << "set xyplane at 0\n";
}

void NeedleAnimation::animate(arma::fvec roa_f_f, arma::fvec euler_angles, 
    std::vector<std::tuple<float, float, float>> beam_points)
{
    create_handle(roa_f_f, euler_angles);
    // create_needle(beam_points);

    // auto xy = beam_points.at(1);
    m_gp << "splot '-' with lines ls 1\n";
    m_gp.send1d(beam_points);
}

void NeedleAnimation::create_needle(std::vector<std::vector<arma::fvec>> beam_points)
{
 
 
}

void NeedleAnimation::create_handle(arma::fvec roa_f_f, arma::fvec euler_angles)
{
    auto cylinder_pol = cylinder(roa_f_f, euler_angles);
    auto str_polygons = polygons_to_strings(cylinder_pol);
    
    int poly_num = 0;

    for (auto str_polyg_i: str_polygons)
    {
        std::string gp_input = "set object " +  std::to_string(poly_num + 1) + 
            " polygon from ";

        int count = 0; 

        for (auto coord: str_polyg_i)
        {
            if (count == str_polyg_i.size() - 1) 
            {
                gp_input += (coord);
            }
            else
            {
                gp_input += (coord + " to ");
            }

            count++; 
        }
        
        gp_input += " fillstyle transparent solid 0.5\n";
        m_gp << gp_input; 
        poly_num++;
    }
        
}


std::vector<arma::fmat> NeedleAnimation::cylinder(arma::fvec center, arma::fvec orientation)
{ 
    
    // Polygons vector
    std::vector<arma::fmat> polygons;
    

    arma::fmat vertices = arma::zeros<arma::fmat>(2 * m_side_count, 3);


    // Cylinder vertices
    for (int i = 0; i < m_side_count; i++)
    {
        float theta = 2 * (M_PIf32 / m_side_count) * i;
        arma::fvec vert_row1 = {-m_height / 2, m_radius * cosf(theta), 
            m_radius * sinf(theta)}; 
        arma::fvec vert_row2 = {m_height / 2, m_radius * cosf(theta), 
            m_radius * sinf(theta)}; 

        // Transform basis 
        arma::fmat rot = EulerRotations().rotation(orientation);
        vert_row1 = center + rot * vert_row1;
        vert_row2 = center + rot * vert_row2;


        vertices.row(i) = vert_row1.t();
        vertices.row(m_side_count + i) = vert_row2.t();
    }

    // Cylinder side faces
    arma::imat side_faces = arma::zeros<arma::imat>(m_side_count, 4);
    for (int i = 0; i < m_side_count - 1; i++)
    {
        arma::ivec side_faces_row = {i + 1, i + 2, m_side_count + i + 2, 
            m_side_count + i + 1}; 
        side_faces.row(i) = side_faces_row.t();
    }
    arma::ivec side_faces_last_row = {m_side_count, 1, m_side_count + 1, 2*m_side_count};
    side_faces.row(m_side_count - 1) = side_faces_last_row.t();

    // Cylinder Bottom faces        
    arma::ivec lower_bottom_face = arma::regspace<arma::ivec>(1, m_side_count);
    arma::ivec upper_bottom_face = arma::regspace<arma::ivec>(m_side_count + 1,  
        2 * m_side_count);

    arma::imat bottom_faces = arma::zeros<arma::imat>(2, m_side_count);
    bottom_faces.row(0) = lower_bottom_face.t();
    bottom_faces.row(1) = upper_bottom_face.t();

    // Vertex processing for plotting 
    for (int i = 0; i < side_faces.n_rows; i++) 
    {
        arma::ivec pol_indices = side_faces.row(i).t();
        arma::fmat polygon_coords = arma::zeros<arma::fmat>(pol_indices.n_rows, 3);

        for (int k = 0; k < pol_indices.n_rows; k++) 
        {
            int index = pol_indices(k);
            polygon_coords.row(k) = vertices.row(index - 1); 
        }
        polygons.push_back(polygon_coords);
    }

    for (int i = 0; i < bottom_faces.n_rows; i++) 
    {
        arma::ivec pol_indices = bottom_faces.row(i).t();
        arma::fmat polygon_coords = arma::zeros<arma::fmat>(pol_indices.n_rows, 3);

        for (int k = 0; k < pol_indices.n_rows; k++) 
        {
            int index = pol_indices(k);
            polygon_coords.row(k) = vertices.row(index - 1); 
        }
        polygons.push_back(polygon_coords);
    }

    return polygons;
}


std::vector<std::vector<std::string>> NeedleAnimation::polygons_to_strings(std::vector<arma::fmat> polygons)
{
    std::vector<std::vector<std::string>> str_polygons;

    for (arma::fmat polyg_coords: polygons)
    {
        std::vector<std::string> coords_vec;

        for (int j = 0; j < polyg_coords.n_rows; j++)
        {
            std::string coords_j = std::to_string(polyg_coords(j, 0)) + "," +
                std::to_string(polyg_coords(j, 1)) + "," + 
                std::to_string(polyg_coords(j, 2));

            coords_vec.push_back(coords_j);
        }

        std::string coords_n = std::to_string(polyg_coords(0, 0)) + "," +
            std::to_string(polyg_coords(0, 1)) + "," + 
            std::to_string(polyg_coords(0, 2));
        coords_vec.push_back(coords_n);
       
        str_polygons.push_back(coords_vec);
    }

    return str_polygons;
}


#endif