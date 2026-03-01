#include <iostream>
#include <string>

#include "Helpers.hpp"
#include "EigenHandler.hpp"
#include "Transformer.hpp"
#include "LinearProblem.hpp"
#include "EnumSolver.hpp"
#include "Solver.hpp"

void print_simplex_solution(const Solver::Solution& res) {
    const auto print_status = [&res]() {
        if (!res.status_message.empty()) {
            std::cout << "Status: " << res.status_message;
        }
        std::cout << std::endl;
    };

    if (!res.is_feasible) {
        std::cout << "Simplex solver failed: task is infeasible" << std::endl;
        print_status();
        return;
    }

    if (res.is_unbounded) {
        std::cout << "Simplex solver failed: task is unbounded" << std::endl;
        print_status();
        return;
    }

    if (!res.is_optimal) {
        std::cout << "Simplex solver failed: result is not optimal" << std::endl;
        print_status();
    }
    std::cout << "Simplex solver solution:" << std::endl;
    std::cout << "\tObjective value: " << res.objective_value << std::endl;
    std::cout << "\tx: ";
    print_vector(res.x);
    std::cout << "\tbasis: ";
    print_vector(res.basis);
    print_status();
}

void print_enum_solution(const EnumSolver::Solution& res) {
    const auto print_status = [&res]() {
        if (!res.status_message.empty()) {
            std::cout << "Status: " << res.status_message;
        }
        std::cout << std::endl;
    };

    if (!res.is_feasible) {
        std::cout << "Enum solver failed: task is infeasible" << std::endl;
        print_status();
        return;
    }

    if (res.is_unbounded) {
        std::cout << "Enum solver failed: task is unbounded" << std::endl;
        print_status();
        return;
    }

    if (!res.is_optimal) {
        std::cout << "Enum solver failed: result is not optimal" << std::endl;
        print_status();
    }
    std::cout << "Enum solver solution:" << std::endl;
    std::cout << "\tObjective value: " << res.objective_value << std::endl;
    std::cout << "\tx: ";
    print_vector(res.x);
    std::cout << "\tbasis: ";
    print_vector(res.basis);
    print_status();
}

int main(int argc, char* argv[]) {
    std::cout << std::fixed << std::setprecision(2);

    LinearProblem original;
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "\nReading problem from file: " << filename << std::endl;
        original = LinearProblem::read_from_file(filename);
    } else {
        std::cout << "\nNo input file specified. Reading from console..." << std::endl;
        original = LinearProblem::read_form_console();
    }

    original.print("Input problem:");

    const auto print_forms = [&original]() {
        LinearProblem general = Transformer::to_general(original);
        general.print("GENERAL FORM (formula 4.1)");
        LinearProblem dual_general = Transformer::to_dual(general);
        dual_general.print("DUAL PROBLEM FOR GENERAL FORM");

        LinearProblem symmetric = Transformer::to_symmetrical(original);
        symmetric.print("SYMMETRIC FORM (formula 4.2)");
        LinearProblem dual_symmetric = Transformer::to_dual(symmetric);
        dual_symmetric.print("DUAL PROBLEM FOR SYMMETRIC FORM");

        LinearProblem canonical = Transformer::to_canonical(original);
        canonical.print("CANONICAL FORM (formula 4.3)");
        LinearProblem dual_canonical = Transformer::to_dual(canonical);
        dual_canonical.print("DUAL PROBLEM FOR CANONICAL FORM");
    };

    // print_forms();

    const auto solve_simplex = [&original]() {
        const auto r = Solver::solve(original, false);
        const auto dual = Transformer::to_dual(original);
        const auto r_dual = Solver::solve(dual, false);
        print_simplex_solution(r);
        print_simplex_solution(r_dual);
        print_val_with_err(r.objective_value, r_dual.objective_value);
    };
    // solve_simplex();

    const auto solve_enum = [&original]() {
        const auto r = EnumSolver::solve(original, false);
        const auto dual = Transformer::to_dual(original);
        const auto r_dual = EnumSolver::solve(dual, false);
        print_enum_solution(r);
        print_enum_solution(r_dual);

        const auto y = std_to_eigen(r_dual.x);
        const auto x = std_to_eigen(r.x);
        const auto A = std_to_eigen(original.constraints());
        const auto c = std_to_eigen(original.target());
// ПРОВЕРКА ОПТИМАЛЬНОСТИ РЕШЕНИЯ ЧЕРЕЗ УСЛОВИЯ КУНА-ТАККЕРА
// 1. Двойственная допустимость: A^T y ≥ c, y ≥ 0
// 2. Дополняющая нежесткость:  (A^T y - c)^T x = 0 (по переменным)
        std::cout << c.transpose() << " - " << y.transpose() * A << std::endl;
        std::cout << c.transpose() - y.transpose() * A << std::endl;
        std::cout << (c.transpose() - y.transpose() * A) << " * " << x << std::endl;
        std::cout << (c.transpose() - y.transpose() * A) * x << std::endl;

        print_val_with_err(r.objective_value, r_dual.objective_value);
    };
    // solve_enum();
    return 0;
}
