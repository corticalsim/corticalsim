#ifndef CSIM_LINALG
#define CSIM_LINALG

#include <Eigen/Dense>

// Symmetric matrix A => eigenvectors in columns of V, corresponding eigenvalues in d
void eigen_decomposition(Eigen::Matrix3d& a, Eigen::Matrix3d& v, Eigen::Vector3d& d);

#endif // CSIM_LINALG
