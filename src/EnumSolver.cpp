#include "EnumSolver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>

#include "Transformer.hpp"
#include "EigenHandler.hpp"
#include "Helpers.hpp"

EnumSolver::Solution EnumSolver::create_infeasible_solution(int n) {
    Solution sol;
    sol.x.resize(n, 0.0);
    sol.objective_value = std::numeric_limits<double>::infinity();
    sol.is_feasible = false;
    sol.is_optimal = false;
    sol.is_unbounded = false;
    sol.is_degenerate = false;
    sol.status_message = "No feasible solution exists";
    return sol;
}

std::vector<std::vector<size_t>> EnumSolver::generate_combinations(int n, int k) {
    std::vector<std::vector<size_t>> result;
    if (k > n || k <= 0)
        return result;

    std::vector<size_t> combination(k);
    std::iota(combination.begin(), combination.end(), 0);

    while (true) {
        result.push_back(combination);

        int i = k - 1;
        while (i >= 0 && combination[i] == static_cast<size_t>(n - k + i)) {
            i--;
        }

        if (i < 0) {
            break; 
        }

        combination[i]++;
        for (int j = i + 1; j < k; ++j) {
            combination[j] = combination[j - 1] + 1;
        }
    }

    return result;
}

bool EnumSolver::is_basis_feasible(const Eigen::VectorXd& xB) {
    return std::ranges::all_of(xB, [](const double d) { return d >= -EPS; });
}

bool EnumSolver::is_basis_degenerate(const Eigen::VectorXd& xB) {
    return std::ranges::any_of(xB, [](const double d) { return std::abs(d) < EPS; });
}

void EnumSolver::print_basis_info(const std::vector<size_t>& basis, const Eigen::VectorXd& xB, double obj_value,
                                  bool degenerate, int solution_index, bool is_new_best, bool verbose) {
    if (!verbose)
        return;

    if (solution_index <= 10 || is_new_best) {
        std::cout << "\nFeasible Basic Solution #" << solution_index << std::endl;

        std::cout << "  Basis indices: ";
        print_vector(basis);

        std::cout << "  Basic variables (x_B): [" << xB << "]" << std::endl;

        if (degenerate) {
            std::cout << "  DEGENERATE SOLUTION (some x_B approx 0)" << std::endl;
        }

        std::cout << "  Objective value: " << obj_value << std::endl;

        if (is_new_best) {
            std::cout << "  >>> NEW BEST SOLUTION <<<" << std::endl;
        }
    } else if (solution_index == 11) {
        std::cout << "\n  ... (showing only first 10 feasible solutions in detail) ..." << std::endl;
    }
}

void EnumSolver::print_final_summary(const Solution& best, size_t total_bases, int feasible_count, int degenerate_count,
                                     int singular_count, bool verbose) {
    if (!verbose) {
        return;
    }

    std::cout << " ENUMERATION COMPLETE: FINAL SUMMARY" << std::endl;

    std::cout << "\nEnumeration statistics:" << std::endl;
    std::cout << "  Total bases examined:       " << total_bases << std::endl;
    std::cout << "  Singular basis matrices:    " << singular_count << std::endl;
    std::cout << "  Feasible basic solutions:   " << feasible_count << std::endl;
    std::cout << "  Degenerate solutions:       " << degenerate_count << std::endl;

    if (!best.is_feasible) {
        std::cout << "\n>>> RESULT: INFEASIBLE PROBLEM <<<" << std::endl;
        std::cout << "  Status: " << best.status_message << std::endl;
        return;
    }

    std::cout << "\nOptimal solution found:" << std::endl;
    std::cout << "  Objective value: " << best.objective_value << std::endl;

    std::cout << "  Solution vector (original variables): ";
    print_vector(best.x);

    std::cout << "  Basis indices: ";
    print_vector(best.basis);

    if (best.is_unbounded) {
        std::cout << "\n  *** WARNING: PROBLEM IS UNBOUNDED ***" << std::endl;
        std::cout << "      Objective can be decreased indefinitely while maintaining feasibility." << std::endl;
    } else if (best.is_optimal) {
        std::cout << "\n  *** SOLUTION IS OPTIMAL ***" << std::endl;
        std::cout << "      All reduced costs are non-negative (minimization problem)." << std::endl;
    }

    if (best.is_degenerate) {
        std::cout << "\n  *** NOTE: SOLUTION IS DEGENERATE ***" << std::endl;
        std::cout << "      At least one basic variable is approximately zero." << std::endl;
    }

    std::cout << "\nStatus: " << best.status_message << std::endl;
}

bool EnumSolver::check_optimality(const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const std::vector<size_t>& basis,
                                  bool& is_unbounded) {
    is_unbounded = false;
    const int m = static_cast<int>(basis.size());
    const int n = static_cast<int>(c.size());

    Eigen::MatrixXd B(m, m);
    Eigen::VectorXd cB(m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            B(i, j) = A(i, static_cast<int>(basis[j]));
        }
        cB(i) = c(static_cast<int>(basis[i]));
    }

    const Eigen::FullPivLU<Eigen::MatrixXd> lu(B);
    if (lu.rank() < m) {
        return false; 
    }

    Eigen::MatrixXd B_inv = lu.inverse();

    Eigen::VectorXd y = B_inv.transpose() * cB;

    for (int j = 0; j < n; ++j) {
        if (std::ranges::find(basis, j) != basis.end()) {
            continue;
        }

        double reduced_cost = c(j) - y.dot(A.col(j));

        if (reduced_cost < -EPS) {
            Eigen::VectorXd direction = B_inv * A.col(j);

            bool can_increase_indefinitely = true;
            for (int i = 0; i < m; ++i) {
                if (direction(i) > EPS) {
                    can_increase_indefinitely = false;
                    break;
                }
            }

            if (can_increase_indefinitely) {
                is_unbounded = true;
                return false; 
            }

            return false;
        }
    }

    return true; 
}

EnumSolver::Solution EnumSolver::evaluate_basis(const Eigen::VectorXd& c, const Eigen::MatrixXd& A,
                                                const Eigen::VectorXd& b, const std::vector<size_t>& basis) {
    const int m = static_cast<int>(basis.size());
    const int n = static_cast<int>(c.size());

    Eigen::MatrixXd B(m, m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            B(i, j) = A(i, static_cast<int>(basis[j]));
        }
    }

    Eigen::FullPivLU<Eigen::MatrixXd> lu(B);
    if (lu.rank() < m || std::abs(lu.determinant()) < EPS) {
        Solution sol;
        sol.is_feasible = false;
        sol.status_message = "Singular basis matrix (rank deficient)";
        return sol;
    }

    Eigen::VectorXd xB;
    try {
        xB = lu.solve(b);
    } catch (...) {
        Solution sol;
        sol.is_feasible = false;
        sol.status_message = "LU decomposition failed";
        return sol;
    }

    if (!is_basis_feasible(xB)) {
        Solution sol;
        sol.is_feasible = false;
        sol.status_message = "Infeasible basis (negative basic variables)";
        return sol;
    }

    Eigen::VectorXd x_full(n);
    x_full.setZero();
    for (int i = 0; i < m; ++i) {
        x_full(static_cast<int>(basis[i])) = xB(i);
    }

    double obj_value = c.dot(x_full);

    Solution sol;
    sol.x.resize(n);
    for (int i = 0; i < n; ++i) {
        sol.x[i] = x_full(i);
    }
    sol.objective_value = obj_value;
    sol.basis = basis;
    sol.is_feasible = true;
    sol.is_degenerate = is_basis_degenerate(xB);
    sol.status_message = sol.is_degenerate ? "Degenerate feasible solution" : "Feasible solution";

    return sol;
}

EnumSolver::Solution EnumSolver::solve(const LinearProblem& lp, bool verbose) {
    if (verbose) {
        lp.print("Enumeration solver. Original problem:");
    }

    LinearProblem canonical = Transformer::to_canonical(lp);

    if (verbose && lp.getForm() != ProblemForm::CANONIC) {
        canonical.print("Problem in canonical form: ");
    }

    int n = static_cast<int>(canonical.num_variables());
    int m = static_cast<int>(canonical.num_constraints());

    Eigen::VectorXd c = std_to_eigen(canonical.target());
    Eigen::MatrixXd A = std_to_eigen(canonical.constraints());
    Eigen::VectorXd b = std_to_eigen(canonical.rhs());

    auto bases = generate_combinations(n, m);
    size_t total_bases = bases.size();

    Solution best = create_infeasible_solution(n);
    best.objective_value = std::numeric_limits<double>::infinity();

    int feasible_count = 0;
    int degenerate_count = 0;
    int singular_count = 0;

    if (verbose) {
        std::cout << "\nEnumerating all basic feasible solutions..." << std::endl;
    }

    for (const auto& basis : bases) {
        Solution sol = evaluate_basis(c, A, b, basis);

        if (!sol.is_feasible) {
            if (sol.status_message.find("Singular") != std::string::npos) {
                singular_count++;
            }
            continue;
        }

        feasible_count++;
        if (sol.is_degenerate) {
            degenerate_count++;
        }

        const bool is_new_best = sol.objective_value < best.objective_value - EPS;

        print_basis_info(basis, Eigen::VectorXd::Map(sol.x.data(), sol.x.size()), sol.objective_value,
                         sol.is_degenerate, feasible_count, is_new_best, verbose && feasible_count <= 10);

        if (is_new_best) {
            best = sol;
            best.is_feasible = true;
        }
    }

    if (!best.is_feasible) {
        best = create_infeasible_solution(lp.num_variables());
        print_final_summary(best, total_bases, feasible_count, degenerate_count, singular_count, verbose);
        return best;
    }

    bool is_unbounded = false;
    best.is_optimal = check_optimality(c, A, best.basis, is_unbounded);
    best.is_unbounded = is_unbounded;

    if (is_unbounded) {
        best.status_message = "Unbounded problem (objective can be decreased indefinitely)";
    } else if (best.is_optimal) {
        best.status_message = "Optimal solution found";
    } else {
        best.status_message = "Feasible solution found (may not be optimal due to numerical issues)";
    }

    best.x = restore_original_solution(lp, best.x);

    if (!lp.type()) {
        best.objective_value = -best.objective_value;
    }

    print_final_summary(best, total_bases, feasible_count, degenerate_count, singular_count, verbose);

    return best;
}
