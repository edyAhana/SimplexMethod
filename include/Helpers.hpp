#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>

#include "LinearProblem.hpp"

std::vector<double> restore_original_solution(const LinearProblem& original_lp,
                                                     const std::vector<double>& canonical_solution);

template <typename T>
void print_vector(const std::vector<T>& v) {
    std::cout << "(";
    if (!v.empty()) {
        for (int i = 0; i < v.size() - 1; i++) {
            std::cout << v[i] << ", ";
        }
        std::cout << v[v.size() - 1];
    }
    std::cout << ")" << std::endl;
}

#endif
