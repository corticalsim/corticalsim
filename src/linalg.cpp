#include "linalg.hpp"

void eigen_decomposition(Eigen::Matrix3d& a, Eigen::Matrix3d& v, Eigen::Vector3d& d)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(a);
    d << solver.eigenvalues();
    v << solver.eigenvectors();
}
