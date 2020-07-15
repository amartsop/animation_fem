#ifndef ARMAEXT_H
#define ARMAEXT_H

#include <iostream>
#include <armadillo>
#include <vector>
#include <fstream>


class ArmaExt
{
    public:
        ArmaExt();

        template <typename T>
        static void sortrows(T *matr, bool clear_reduntant=false);

    private:
        template <typename T>
        static bool compare_head(const arma::Mat<typename T::elem_type>& lhs, const arma::Mat<typename T::elem_type>& rhs);
};


/**
 * Sorts in ascending order an Armadillo matrix of type T according to it's first row 
 * @param mat Armadillo matrix of type T
 * @return Sorted Armadillo matrix of type T
*/
template <typename T>
void ArmaExt::sortrows(T *matr, bool clear_reduntant)
{  
    std::vector< arma::Mat<typename T::elem_type> > vec;


    for (int64_t i = 0; i < matr->n_rows; ++i)
    {
        vec.push_back(matr->row(i));
    }

    std::sort(vec.begin(), vec.end(), &compare_head<T>);

    for (int64_t i = 0; i < matr->n_rows; ++i)
    {
            matr->row(i) = vec[i];
    }

    if(clear_reduntant)
    {
        int64_t matu_rows = 0;
        for (int64_t i = 0; i < matr->n_rows - 1; i++)
        {
            if ((*matr)(i + 1, 0) - (*matr)(i, 0) != 0) { matu_rows++ ;}
        }

        arma::Mat<typename T::elem_type> mat_u(matu_rows + 1, (*matr).n_cols); 
        mat_u.zeros();        
        int64_t counter = 0;
        
        for (int64_t i = 0; i < matr->n_rows - 1; i++)
        {
            if ((*matr)(i + 1, 0) - (*matr)(i, 0) != 0) 
            { 
                mat_u.row(counter) = matr->row(i);
                counter++;
            }
        }
        
        mat_u.row(mat_u.n_rows - 1) = matr->row(matr->n_rows - 1);
        (*matr) = mat_u;
    }
    
}

template <typename T>
bool ArmaExt::compare_head(const arma::Mat<typename T::elem_type>& lhs, const arma::Mat<typename T::elem_type>& rhs)
{
    return lhs(0) < rhs(0);
}


#endif