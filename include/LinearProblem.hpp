#ifndef LINEAR_PROBLEM_H
#define LINEAR_PROBLEM_H

#include <vector>
#include <string>
#include <filesystem>

using std::vector;
using std::string;
using std::filesystem::path;

enum class ProblemForm {
    GENERAL,
    CANONIC, 
    SYMMETRIC,
};

class LinearProblem {
private:
    bool problem_type_;
    vector<double> target_;
    vector<vector<double>> constraints_;
    vector<string> constraint_type_;
    vector<double> rhs_;
    vector<string> var_constraints_;

    static string form_to_string(ProblemForm form);
public:
    LinearProblem() = default;
    LinearProblem(bool problem_type,
                  vector<double> target,
                  vector<vector<double>> constraints,
                  vector<string> constraint_type,
                  vector<double> rhs,
                  vector<string> var_constraints);

    static LinearProblem read_from_file(const path& file);
    static LinearProblem read_form_console();

    // getters
    
    auto num_constraints() { return constraints_.size(); }
    auto num_variables() { return target_.size(); }

    auto& type() { return problem_type_; }
    auto& target() { return target_; }
    auto& constraints() { return constraints_; }
    auto& constraint_type() { return constraint_type_; }
    auto& rhs() { return rhs_; }
    auto& var_constaints() { return var_constraints_; }

    const auto& type() const { return problem_type_; }
    const auto& target() const { return target_; }
    const auto& constraints() const { return constraints_; }
    const auto& constraint_type() const { return constraint_type_; }
    const auto& rhs() const { return rhs_; }
    const auto& var_constaints() const { return var_constraints_; }

    auto getM1() const;
    auto getM2() const;
    auto getN1() const;
    auto getN2() const;

    auto getForm() const;

    void print(const string& title) const;
};



#endif
