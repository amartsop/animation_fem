#ifndef EULER_ROTATIONS_H
#define EULER_ROTATIONS_H

#include <iostream>
#include <armadillo>

class EulerRotations
{
public:
    EulerRotations();
    static arma::fmat basic_rotation_x(float x);
    static arma::fmat basic_rotation_y(float x);
    static arma::fmat basic_rotation_z(float x);
    static arma::fmat rotation(float phi, float theta, float psi);
    static arma::fmat rotation(arma::fvec euler_angles);

private:

};


EulerRotations::EulerRotations()
{
}

arma::fmat EulerRotations::basic_rotation_x(float x)
{
    arma::fmat rot = {
        {1.0f, 0.0f, 0.0f}, 
        {0.0f, cosf(x), -sinf(x)}, 
        {0.0f, sinf(x), cosf(x)}};
    return rot;
}

arma::fmat EulerRotations::basic_rotation_y(float x)
{
    arma::fmat rot = {
        {cosf(x), 0.0f, sinf(x)}, 
        {0.0f, 1.0f, 0.0f}, 
        {-sinf(x), 0.0f, cosf(x)}};
    return rot;
}

arma::fmat EulerRotations::basic_rotation_z(float x)
{
    arma::fmat rot = {
        {cosf(x), -sinf(x), 0.0f}, 
        {sinf(x), cosf(x), 0.0f}, 
        {0.0f, 0.0f, 1.0f}};
    return rot;
}

arma::fmat EulerRotations::rotation(float phi, float theta, float psi)
{

    arma::fmat rotx =  basic_rotation_x(phi);
    arma::fmat roty = basic_rotation_y(theta);
    arma::fmat rotz = basic_rotation_z(psi);

    arma::fmat rot = (rotz * roty * rotx);
    return rot;
}

arma::fmat EulerRotations::rotation(arma::fvec euler_angles)
{
    return rotation(euler_angles(0), euler_angles(1), euler_angles(2));
}

#endif