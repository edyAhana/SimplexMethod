#include <iostream>
#include <stdexcept>

#include "Transformer.hpp"

LinearProblem Transformer::to_general(const LinearProblem& lp) {
    // if (lp.getForm() == ProblemForm::GENERAL) {
    //     return lp;
    // }
    
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
    bool problem_type = lp.type();
    auto target = lp.target();
    auto constraints = lp.constraints();
    auto constraint_type = lp.constraint_type();
    auto rhs = lp.rhs();
    auto var_constraints = lp.var_constaints();

    auto var_num = target.size();
    auto constraint_num = constraints.size();

    auto N2 = lp.getN2();

    if(!problem_type) {
        for(auto coef : target) {
            coef *= -1.;
        }
        problem_type = true;
    }

    for(std::size_t i = 0; i < constraint_num; ++i) {
        if(constraint_type[i] == "<=") {
            for(std::size_t j = 0; j < var_num; ++j) {
                constraints[i][j] *= -1;
            }
            rhs[i] *= -1;
            constraint_type[i] = ">=";
        }
    }

    
    for(std::size_t i = 0; i < N2.size(); ++i) {
        if(var_constraints[N2[i]] == "<=") {
            for(std::size_t j = 0; j < constraint_num; ++j) {
                constraints[j][N2[i]] *= -1;
            }
            constraint_type[N2[i]] = ">=";
        }
    }

    return LinearProblem(problem_type, target, constraints,
                        constraint_type, rhs, var_constraints);
}

LinearProblem Transformer::to_canonical_form(const LinearProblem& lp) {
    auto& target = lp.target();
    auto& constraints = lp.constraints();
    auto& constraint_type = lp.constraint_type();
    auto rhs = lp.rhs();
    auto& var_constraints = lp.var_constaints();
    
    size_t num_vars = target.size();
    size_t num_cons = constraints.size();
    
    auto M1 = lp.getM1();
    auto M2 = lp.getM2();
    auto N1 = lp.getN1();
    auto N2 = lp.getN2();

    auto new_num_var = N1.size() + N2.size()*2 + M1.size();

    vector<double> new_target(new_num_var);
    vector<vector<double>> new_constraints(num_cons, vector<double>(new_num_var));
    vector<double> new_rhs(num_cons);

    std::size_t var_cnt = 0;

    for(std::size_t i = 0; i < num_vars; ++i) {
        if(var_constraints[i] == ">=0") {
            new_target[var_cnt++] = target[i];
        } else {
            new_target[var_cnt++] = target[i];
            new_target[var_cnt++] = -1 * target[i];
        }
    }

    var_cnt = 0;

    for(std::size_t i = 0; i < num_vars; ++i) {
        if(var_constraints[i] == ">=0") {
            for(std::size_t j = 0; j < num_cons; ++j) {
                new_constraints[j][var_cnt] = constraints[j][i];
            }
            ++var_cnt;
        } else {
            for(std::size_t j = 0; j < num_cons; ++j) {
                new_constraints[j][var_cnt] = constraints[j][i];
                new_constraints[j][var_cnt+1] = -1*constraints[j][i];
            }
            var_cnt += 2;
        }
    }
    
    int slack_index = N1.size() + N2.size()*2;

    for (size_t i = 0; i < num_cons; ++i) {
        if (constraint_type[i] == ">=") {
            new_constraints[i][slack_index] = -1.0;
            slack_index++;
        } 
    }
    
    return LinearProblem(lp.type(), new_target, new_constraints, 
                        vector<string>(num_cons, "="), rhs, vector<string>(new_num_var, ">=0"));
}

LinearProblem Transformer::to_symmetrical_form(const LinearProblem& lp) {
    auto& target = lp.target();
    auto& constraints = lp.constraints();
    auto& constraint_type = lp.constraint_type();
    auto& rhs = lp.rhs();
    auto& var_constraints = lp.var_constaints();
    
    size_t num_vars = target.size();
    size_t num_cons = constraints.size();
    
    auto M1 = lp.getM1();
    auto M2 = lp.getM2();
    auto N1 = lp.getN1();
    auto N2 = lp.getN2();

    auto new_num_var = N1.size() + N2.size()*2;
    auto new_num_cons = M1.size() + M2.size()*2;

    vector<double> new_target(new_num_var);
    vector<vector<double>> new_constraints(new_num_cons, vector<double>(new_num_var));
    vector<double> new_rhs(new_num_cons);

    std::size_t var_cnt = 0;

    for(std::size_t i = 0; i < num_vars; ++i) {
        if(var_constraints[i] == ">=0") {
            new_target[var_cnt++] = target[i];
        } else {
            new_target[var_cnt++] = target[i];
            new_target[var_cnt++] = -1 * target[i];
        }
    }

    var_cnt = 0;
    for(std::size_t i = 0; i < num_vars; ++i) {
        std::size_t cons_cnt = 0;
        for(std::size_t j = 0; j < num_cons; ++j) {
            if(constraint_type[j] == ">=") {
                if(var_constraints[i] == ">=0") {
                    new_constraints[cons_cnt][var_cnt] = constraints[j][i];
                    ++cons_cnt;
                } else {
                    new_constraints[cons_cnt][var_cnt] = constraints[j][i];
                    new_constraints[cons_cnt][var_cnt+1] = -1 * constraints[j][i];
                    ++cons_cnt;
                }
            } else {
                if(var_constraints[i] == ">=0") {
                    new_constraints[cons_cnt][var_cnt] = constraints[j][i];
                    new_constraints[cons_cnt+1][var_cnt] = -1 *constraints[j][i];
                    cons_cnt+=2;
                } else {
                    new_constraints[cons_cnt][var_cnt] = constraints[j][i];
                    new_constraints[cons_cnt][var_cnt+1] = -1 * constraints[j][i];
                    new_constraints[cons_cnt+1][var_cnt] = -1 * constraints[j][i];
                    new_constraints[cons_cnt+1][var_cnt+1] =  constraints[j][i];
                    cons_cnt+=2;
                }
            }
        }
        var_cnt = var_constraints[i] == ">=0" ? var_cnt+1: var_cnt+2;
    }

    std::size_t cons_cnt = 0;
    for(std::size_t i = 0; i < num_cons; ++i) {
        if(constraint_type[i] == ">=") {
            new_rhs[cons_cnt] = rhs[i];
            ++cons_cnt;
        } else {
            new_rhs[cons_cnt] = rhs[i];
            new_rhs[cons_cnt+1] = -1*rhs[i];
            cons_cnt += 2;
        }
    }

    return LinearProblem(lp.type(), new_target, new_constraints, 
                        vector<string>(new_num_cons, ">="), new_rhs, vector<string>(new_num_var, ">=0"));
    
}


