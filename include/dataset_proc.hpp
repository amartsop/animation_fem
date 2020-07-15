#ifndef DATASET_PROCESSING_H
#define DATASET_PROCESSING_H

#include <armadillo>

class DatasetProcessing
{
    public:
        DatasetProcessing();

       static void resampling(arma::fmat *matr, float ts);

    private:
};

/**
 * Resamples a sorted Armadillo matrix (ascenting) of type fmat in  intervals of ts 
 * (using linear interpolation).
 * @param matr Armadillo matrix of type fmat
 * @return Rasampled Armadillo matrix of type fmat
*/
void DatasetProcessing::resampling(arma::fmat *matr, float ts)
{
    arma::fvec t = (*matr).col(0);
    arma::fvec tu = arma::regspace<arma::fvec>(t.front(), ts, t.back());
    arma::fmat mat_u(tu.n_rows, (*matr).n_cols); mat_u.zeros();
    mat_u(0, arma::span::all) = (*matr)(0, arma::span::all);

    for (int64_t i = 1; i < mat_u.n_rows; i++)
    {   
        arma::uvec q1 =  arma::find(t < tu(i));
        arma::uvec q2 = arma::find(t > tu(i));

        // Indices that tu(i) belongs in
        int64_t index1 = q1.back(); int64_t index2 = q2.front();

        // Times that tu(i) belongs in
        float t1 = t(index1); float t2 = t(index2);

        // Corresponding values for t1 and t2
        arma::fvec a1 = (*matr)(index1, arma::span(1, matr->n_cols - 1));
        arma::fvec a2 = (*matr)(index2, arma::span(1, matr->n_cols - 1));

        // Interpolation
        arma::fvec a = a2 - ((a2 - a1) / (t2 - t1)) * (t2 - tu(i));

        mat_u(i, 0) = tu(i);
        mat_u(i, arma::span(1, mat_u.n_cols - 1)) = a;
    }
    
    (*matr) = mat_u;
}

#endif