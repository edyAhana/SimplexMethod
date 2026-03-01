#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

#include "LinearProblem.hpp"

class EnumSolver {
public:
    struct Solution {
        std::vector<double> x;
        double objective_value;
        std::vector<size_t> basis;
        bool is_feasible;
        bool is_optimal;
        bool is_unbounded;
        bool is_degenerate;
        std::string status_message;
    };

    static Solution solve(const LinearProblem& lp, bool verbose = false);

private:
    static constexpr double EPS = 1e-9;

    static Solution create_infeasible_solution(int n);
    static std::vector<std::vector<size_t>> generate_combinations(int n, int k);
    static bool is_basis_feasible(const Eigen::VectorXd& xB);
    static bool is_basis_degenerate(const Eigen::VectorXd& xB);

    static Solution evaluate_basis(const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                   const std::vector<size_t>& basis);
    static void print_basis_info(const std::vector<size_t>& basis, const Eigen::VectorXd& xB, double obj_value,
                                 bool degenerate, int solution_index, bool is_new_best, bool verbose);
    static void print_final_summary(const Solution& best, size_t total_bases, int feasible_count, int degenerate_count,
                                    int singular_count, bool verbose);
    static bool check_optimality(const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const std::vector<size_t>& basis,
                                 bool& is_unbounded);
};
