#include "linalg.hpp"

void eigen_decomposition(Eigen::Matrix3d& A, Eigen::Matrix3d& V, Eigen::Vector3d& d)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(A);
    d << solver.eigenvalues();
    V << solver.eigenvectors();
}
