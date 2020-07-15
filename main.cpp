#include <iostream>
#include <armadillo>
#include <vector>
#include <chrono>
#include "needle_animation.hpp"
#include "euler_rotations.hpp"
#include "needle_geometry.hpp"

arma::fmat elastic_coordinates(float t, int ne);
arma::fmat shape_function(float x2j, float le);
arma::fmat locator_matrix(int element_num, int elastic_coordinates_num);
arma::fvec rigid_body_position(float t);
arma::fvec rigid_body_orientation(float t);
std::vector<std::tuple<float, float, float>> beam_coordinates(float t,
    arma::fvec roa_f_f, arma::fvec euler_angles);

int main(int argc, char **argv) {
    
    auto start = std::chrono::steady_clock::now();
    
    float t = 0;

    NeedleAnimation na;
    
    while(1)
    {
   
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<float> elapsed_seconds = end-start;
    
    // Update time
    // Real time
    float t = t + elapsed_seconds.count();
    
    // Update body position 
    auto roa_f_f = rigid_body_position(t);
    auto euler_angles = rigid_body_orientation(t);
    auto beam_points = beam_coordinates(t, roa_f_f, euler_angles);

    // Animate
    na.animate(roa_f_f, euler_angles, beam_points);

    // Update time
    start = end;

    }

}

arma::fvec rigid_body_position(float t)
{
    float fs = 1;
    float sint = sinf(2 * M_PIf32 * t);
    float cost = cosf(2 * M_PIf32 * t);

    arma::fvec roa_f_f = {0.05f * sint, 0.1f * cost, 0.04f * sint};

    return roa_f_f;

}

arma::fvec rigid_body_orientation(float t)
{
    float fs = 1;
    float sint = sinf(2 * M_PIf32 * t);
    float cost = cosf(2 * M_PIf32 * t);

    arma::fvec euler_angles = {0.01f * cost, 0.3f * sint, 0.5f * cost};

    return euler_angles;
}

std::vector<std::tuple<float, float, float>> beam_coordinates(float t,
    arma::fvec roa_f_f, arma::fvec euler_angles)
{

    // Number of elements
    int ne = 10; 
    
    // Number of dofs (elastic coordinates)
    int n = 2 * ne + 2;
    
    // Rigid body coordinates
    arma::fmat rot_f1_f = EulerRotations::rotation(euler_angles);
    float lx = NeedleGeometry::get_needle_offset_x();
    float lz = NeedleGeometry::get_needle_offset_z();
    arma::fvec rab_f1_f1 = {lx, 0, -lz};
    arma::fvec rob_f_f = roa_f_f + rot_f1_f * rab_f1_f1;

    // Elastic coordinates
    arma::fvec qe = elastic_coordinates(t, ne);

    // Beam length
    float l = NeedleGeometry::get_needle_length();

    // Element length
    float le = l / ne;
    
    // Discretize element
    int npel = 10; // Number of point in each element
    arma::fvec x2j_vec = arma::linspace<arma::fvec> (0, le, npel);

    // Positions of the flexible beam in time ti
    std::vector<std::tuple<float, float, float>> beam_points;

    // Shape functions
    for (int j = 1; j <= ne; j++)
    {
        arma::fmat lj_mat = locator_matrix(j, n);

        for (int i = 0; i < npel; i++) 
        {
            arma::fmat psi_mat = shape_function(x2j_vec(i), le);

            // Initial position of pj wrt to f2 origin
            arma::fvec rbp0j_f2_f2 = { (float) (j - 1) * le + x2j_vec(i), 0.0f, 0.0f};

            // Deformation approximation
            arma::fvec rp0pj_f2_f2 = psi_mat * lj_mat * qe;

            // Final position of pj wrt to f2 origin
            arma::fvec rbpj_f2_f2 = rbp0j_f2_f2 + rp0pj_f2_f2;

            // Final position of pj wrt to F origin 
            arma::fvec ropj_f_f = rob_f_f + rot_f1_f *  rbpj_f2_f2;

            // Final position of pj wrt to F origin tuple
            std::tuple ropj_f_f_tuple = std::make_tuple(ropj_f_f(0), ropj_f_f(1),
                ropj_f_f(2));

            beam_points.push_back(ropj_f_f_tuple);
        }
    }
    return beam_points;
}


arma::fmat elastic_coordinates(float t, int ne)
{
    arma::fvec qe(2 * ne + 2); qe.fill(1);
    arma::fvec amplitude(2 * ne + 2); amplitude.fill(0.01);
    
    arma::fvec phase(2 * ne + 2); phase.fill(0);
    arma::fvec frequency = arma::linspace<arma::fvec>(0, 100, 2 * ne + 2);

    for(int i = 0; i < (2 * ne + 2); i++) 
    {
        qe(i) = amplitude(i) * sinf(2 * M_PIf32 * frequency(i) * t + phase(i));
    }

    return qe;
}


arma::fmat shape_function(float x2j, float le)
{
    float ksi = x2j / le;

    arma::fmat psi = arma::zeros<arma::fmat>(3, 4);
    psi(2, 0) = 1 - 3 * powf(ksi, 2) + 2 * powf(ksi, 3);
    psi(2, 1) = x2j - 2 * le * powf(ksi, 2) + le * powf(ksi, 3);
    psi(2, 2) = 3 * powf(ksi, 2) - 2 * powf(ksi, 3);
    psi(2, 3) = -le * powf(ksi, 2) + le * powf(ksi, 3);

    return psi;
}

arma::fmat locator_matrix(int element_num, int elastic_coordinates_num)
{
    arma::fmat l_mat = arma::zeros<arma::fmat>(4, elastic_coordinates_num);

    // 2D Beam
    int index1 = (2 * element_num - 1) - 1; 
    l_mat(0, index1) = 1;
   
    int index2 = (2 * element_num) - 1; 
    l_mat(1, index2) = 1;
   
    int index3 = (2 * element_num + 1) - 1; 
    l_mat(2, index3) = 1;
    
    
    int index4 = (2 * element_num + 2) - 1; 
    l_mat(3, index4) = 1;

    return l_mat;
}