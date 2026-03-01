#ifndef SOLVER_H
#define SOLVER_H

#include <optional>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "LinearProblem.hpp"

class Solver {
private:
    static constexpr std::size_t MAX_ITER = 1000;
    static constexpr double EPS = 10e-9;

    struct PhaseResult {
            bool optimal;
            bool unbounded;
            int iterations;
            Eigen::MatrixXd B_inv;
            Eigen::VectorXd x_B;
            std::vector<size_t> basis;
            double objective_value;
    };
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
private:
    static std::optional<PhaseResult> run_simplex_phase(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const vector<double> c,
                                  std::vector<std::size_t> basis, Eigen::MatrixXd B_inv, Eigen::VectorXd x_B,
                                  std::size_t n_active_vars, string phase_name, bool verbose);

    static bool cleanup_artificial_basis(const Eigen::MatrixXd& A_orig, Eigen::MatrixXd& B_inv,
                                         std::vector<size_t>& basis, Eigen::VectorXd& x_B, size_t n_orig, bool verbose);

    static Solution create_infeasible_solution(int n);

    static bool is_basis_feasible(const Eigen::VectorXd& xB, double tolerance);

    static bool is_basis_degenerate(const Eigen::VectorXd& xB, double tolerance);

    static void print_basis_change(const std::string& phase_name, int iteration, size_t entering_var,
                                   size_t leaving_var, size_t leaving_pos, double reduced_cost, double theta,
                                   bool degenerate, const std::vector<size_t>& new_basis, const Eigen::VectorXd& new_xB,
                                   double objective_value, bool verbose);

    static void print_final_summary(const Solution& sol, int phase1_iterations, int phase2_iterations, bool verbose);

    static Solution handle_no_constraints(size_t n, const std::vector<double>& c, bool verbose);

    static Solution handle_no_variables(size_t m, const std::vector<double>& b, bool verbose);

    static bool check_optimality(const Eigen::MatrixXd& A, const Eigen::VectorXd& c, const std::vector<size_t>& basis,
                                 const Eigen::MatrixXd& B_inv, bool& is_unbounded);
public:
    static Solution solve(const LinearProblem& lp, bool verbose);
};


#endif
