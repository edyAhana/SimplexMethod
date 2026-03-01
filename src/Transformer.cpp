#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "Transformer.hpp"

LinearProblem Transformer::to_general(const LinearProblem& lp) {
    if (lp.getForm() == ProblemForm::GENERAL) {
        return lp;
    }
    
    return to_general_form(lp);
}

LinearProblem Transformer::to_canonical(const LinearProblem& lp) {
    if (lp.getForm() == ProblemForm::CANONIC) {
        return lp;
    }
    
    LinearProblem general_form = to_general(lp);
    return to_canonical_form(general_form);
}

LinearProblem Transformer::to_symmetrical(const LinearProblem& lp) {
    if (lp.getForm() == ProblemForm::SYMMETRIC) {
        return lp;
    }
    
    LinearProblem general_form = to_general(lp);
    return to_symmetrical_form(general_form);
}

LinearProblem Transformer::to_general_form(const LinearProblem& lp) {
    LinearProblem general = lp;
    for (size_t i = 0; i < lp.num_variables(); ++i) {
        if (lp.var_constaints()[i] == "<=0") {
            general.target()[i] = -lp.target()[i];
            for (size_t j = 0; j < lp.num_constraints(); ++j) {
                general.constraints()[i][j] = -lp.constraints()[i][j];
            }

            general.var_constaints()[i] = ">=0";
        }
    }

    // max -> min
    if (!lp.type()) {
        for (double& coeff : general.target())
            coeff = -coeff;
        general.type() = true;
    }

    std::vector<std::vector<double>> new_constraints;
    std::vector<std::string> new_relations;
    std::vector<double> new_rhs;

    for (size_t i = 0; i < lp.num_constraints(); ++i) {
        if (lp.constraint_type()[i] == "<=") {
            // a[i]x <= b[i] -> a[i]x >= b[i]
            std::vector<double> neg_row(lp.num_variables());
            for (size_t j = 0; j < lp.num_variables(); ++j) {
                neg_row[j] = -lp.constraints()[i][j];
            }
            new_constraints.emplace_back(neg_row);
            new_relations.emplace_back(">=");
            new_rhs.emplace_back(-lp.rhs()[i]);
        } else if (lp.constraint_type()[i] == ">=" || lp.constraint_type()[i] == "=") {
            new_constraints.emplace_back(lp.constraints()[i]);
            new_relations.emplace_back(lp.constraint_type()[i]);
            new_rhs.emplace_back(lp.rhs()[i]);
        } else {
            throw std::runtime_error("Invalid constraint type: " + lp.constraint_type()[i]);
        }
    }

    general.constraints() = std::move(new_constraints);
    general.constraint_type() = std::move(new_relations);
    general.rhs() = std::move(new_rhs);
    general.prettierCoeffs();

    return general;
}

LinearProblem Transformer::to_canonical_form(const LinearProblem& lp) {
    LinearProblem general = to_general_form(lp);

    auto& constraints = general.constraints();
    auto& rhs = general.rhs();

    const auto M1 = general.getM1();
    const auto M2 = general.getM2();
    const auto N1 = general.getN1();
    const auto N2 = general.getN2();

    size_t new_vars = N2.size() + M1.size();

    // c = (c[N1], c[N2], -c[N2], 0[M1])
    for (const auto& i : N2) {
        general.target().emplace_back(-general.target()[i]);
    }
    for (int i = 0; i < M1.size(); i++) {
        general.target().emplace_back(0);
    }

    // A = (A[M1,N1], A[M1,N2], -A[M1,N2], -E[M1,M1]
    //     (A[M2,N1], A[M2,N2], -A[M2,N2],  0[M2,M1]
    for (size_t j = 0; j < M1.size(); j++) {
        const size_t y = M1[j];
        for (const auto& i : N2) {
            constraints[y].emplace_back(-general.constraints()[y][i]);
        }
        for (size_t i = 0; i < M1.size(); i++) {
            if (i == j) {
                constraints[y].emplace_back(-1);
            } else {
                constraints[y].emplace_back(0);
            }
        }
    }
    for (const auto& y : M2) {
        for (const auto& i : N2) {
            constraints[y].emplace_back(-general.constraints()[y][i]);
        }
        for (size_t i = 0; i < M1.size(); i++) {
            constraints[y].emplace_back(0);
        }
    }

    std::ranges::fill(general.var_constaints(), ">=0");
    for (size_t i = 0; i < new_vars; i++) {
        general.var_constaints().emplace_back(">=0");
    }
    // now new constraints
    std::ranges::fill(general.constraint_type(), "=");

    // make rhs >= 0
    for (size_t i = 0; i < general.num_constraints(); ++i) {
        if (rhs[i] < 0) {
            for (size_t j = 0; j < general.num_variables(); ++j) {
                constraints[i][j] = -constraints[i][j];
            }
            rhs[i] = -rhs[i];
        }
    }

    general.prettierCoeffs();
    return general;
}

LinearProblem Transformer::to_symmetrical_form(const LinearProblem& lp) {
    LinearProblem general = to_general_form(lp);

    const auto M1 = general.getM1();
    const auto M2 = general.getM2();
    const auto N1 = general.getN1();
    const auto N2 = general.getN2();

    size_t original_vars = lp.num_variables();
    size_t new_vars = N2.size();
    size_t total_vars = original_vars + new_vars;

    size_t original_constraints = general.num_constraints();
    size_t new_constraints = M2.size();
    size_t total_constraints = original_constraints + new_constraints;

    // c = (c[N1], c[N2], -c[N2])
    for (const auto& i : N2) {
        general.target().emplace_back(-general.target()[i]);
    }

    // b = (b[M1], b[M2], -b[M2])
    for (const auto& i : M2) {
        general.rhs().emplace_back(-general.rhs()[i]);
    }


    auto& constraints = general.constraints();

    //     ( A[M1][N1],  A[M1][N2], -A[M1,N2])
    // A = ( A[M2][N1],  A[M2][N2], -A[M2,N2])
    //     (-A[M2][N1], -A[M2][N2],  A[M2,N2])
    for (size_t y = 0; y < constraints.size(); y++) {
        for (const auto& i : N2) {
            constraints[y].emplace_back(-constraints[y][i]);
        }
    }
    for (size_t y = 0; y < M2.size(); y++) {
        constraints.emplace_back();
        constraints[y + original_constraints].reserve(total_vars);
        for (size_t x = 0; x < original_vars; x++) {
            constraints[y + original_constraints].emplace_back(-constraints[M2[y]][x]);
        }
        for (const auto& i : N2) {
            constraints[y + original_constraints].emplace_back(constraints[M2[y]][i]);
        }
    }

    std::ranges::fill(general.constraint_type(), ">=");
    for (size_t i = 0; i < new_constraints; i++) {
        general.constraint_type().emplace_back(">=");
    }

    std::ranges::fill(general.var_constaints(), ">=0");
    for (size_t i = 0; i < new_vars; i++) {
        general.var_constaints().emplace_back(">=0");
    }

    return general;
}

LinearProblem Transformer::to_dual(const LinearProblem& problem) {
    auto lp = to_general(problem);
    auto& target = lp.target();
    auto& constraints = lp.constraints();
    auto& constraint_type = lp.constraint_type();
    auto& rhs = lp.rhs();
    auto& var_constraints = lp.var_constaints();

    auto new_target = rhs;
    auto new_rhs = target;

    auto num_vars = lp.num_variables();
    auto num_cons = lp.num_constraints();

    vector<vector<double>> new_constraints(num_vars, vector<double>(num_cons));
    vector<string> new_constraint_type(num_vars);
    vector<string> new_var_constraints(num_cons);
    
    bool new_type = !lp.type();



    for(std::size_t i = 0; i < num_cons; ++i) {
        for(std::size_t j = 0; j < num_vars; ++j) {
            new_constraints[j][i] = constraints[i][j];
        }
    }

    const auto M1 = lp.getM1();
    const auto M2 = lp.getM2();
    const auto N1 = lp.getN1();
    const auto N2 = lp.getN2();

    for(std::size_t i = 0; i < N1.size(); ++i) {
        new_constraint_type[N1[i]] = "<=";
    }

    for(std::size_t i = 0; i < N2.size(); ++i) {
        new_constraint_type[N2[i]] = "=";
    }

    for(std::size_t i = 0; i < M1.size(); ++i) {
        new_var_constraints[M1[i]] = ">=0";
    }

    for(std::size_t i = 0; i < M2.size(); ++i) {
        new_var_constraints[M2[i]] = "free";
    }

    return LinearProblem(new_type, new_target, new_constraints, new_constraint_type, new_rhs, new_var_constraints);
}




