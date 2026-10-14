#include "linalg.hpp"
#include "randhub.hpp"
#include "eig3.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <limits>
#include <random>

template <typename T = double[]>
void print_mat(T mat, std::string header)
{
    std::cout << header << std::endl;
    for (size_t r = 0; r < 3; ++r)
    {
        for (size_t c = 0; c < 3; ++c)
        {
            if constexpr (std::is_same_v<T, Eigen::Matrix3d>)
            {

                std::cout << mat(r, c) << '\t';
            }
            else
            {

                std::cout << mat[r][c] << '\t';
            }
        }
        std::cout << std::endl;
    }
}

template <typename T = double[]>
void print_vec(T vec, std::string header)
{
    std::cout << header << std::endl;
    for (size_t elem = 0; elem < 3; ++elem)
    {
        std::cout << vec[elem] << std::endl;
    }
}

int main(int argc, char* argv[])
{

    RandHub<> rhub = RandHub();

    double qxx = rhub.rand(-1, 1);
    double qxy = rhub.rand(-1, 1);
    double qxz = rhub.rand(-1, 1);
    double qyy = rhub.rand(-1, 1);
    double qyz = rhub.rand(-1, 1);
    double qzz = rhub.rand(-1, 1);

    double matrix[3][3] = { { qxx, qxy, qxz }, { qxy, qyy, qyz }, { qxz, qyz, qzz } };
    double evecs[3][3];
    double evals[3];

    std::cout << std::setprecision(3);

    print_mat<>(matrix, "\n==[ double[] matrix ]==");

    eigen_decomposition_manual(matrix, evecs, evals);

    print_mat<>(evecs, "\n==[ double[] eigenvectors ]==");
    print_vec<>(evals, "\n==[ double[] eigenvalues ]==");

    // ==========================================

    Eigen::Matrix3d eigen_matrix{ { qxx, qxy, qxz }, { qxy, qyy, qyz }, { qxz, qyz, qzz } };
    Eigen::Matrix3d eigen_evecs{ Eigen::Matrix3d::Zero() };
    Eigen::Vector3d eigen_evals{ Eigen::Vector3d::Zero() };

    print_mat<>(eigen_matrix, "\n==[ Eigen matrix ]==");

    eigen_decomposition(eigen_matrix, eigen_evecs, eigen_evals);

    print_mat<>(eigen_evecs, "\n==[ Eigen eigenvectors ]==");
    print_vec<>(eigen_evals, "\n==[ Eigen eigenvalues ]==");

    return 0;
}
