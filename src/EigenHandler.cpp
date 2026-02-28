#include "EigenHandler.hpp"

Eigen::VectorXd std_to_eigen(const std::vector<double>& v) {
    Eigen::VectorXd result;
    result.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i];
    }
    return result;
}

Eigen::VectorXd std_to_eigen(const std::vector<double>& v, const std::function<double(double)>& transform) {
    Eigen::VectorXd result;
    result.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = transform(v[i]);
    }
    return result;
}

Eigen::MatrixXd std_to_eigen(const std::vector<std::vector<double>>& m,
                                    const std::function<double(double)>& transform) {
    Eigen::MatrixXd result;
    result.resize(m.size(), m[0].size());
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[i].size(); ++j) {
            result(i, j) = transform(m[i][j]);
        }
    }
    return result;
}

Eigen::MatrixXd std_to_eigen(const std::vector<std::vector<double>>& m) {
    Eigen::MatrixXd result;
    result.resize(m.size(), m[0].size());
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[i].size(); ++j) {
            result(i, j) = m[i][j];
        }
    }
    return result;
}
