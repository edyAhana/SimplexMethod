#ifndef EIGEN_HANDLER_H
#define EIGEN_HANDLER_H

#include <Eigen/Dense>
#include <cassert>
#include <random>
#include <vector>

Eigen::VectorXd std_to_eigen(const std::vector<double>& v);

Eigen::VectorXd std_to_eigen(const std::vector<double>& v,
                                    const std::function<double(double)>& transform);

Eigen::MatrixXd std_to_eigen(const std::vector<std::vector<double>>& m);

Eigen::MatrixXd std_to_eigen(const std::vector<std::vector<double>>& m,
                                    const std::function<double(double)>& transform);



#endif
