#ifndef LINEAR_PROBLEM_H
#define LINEAR_PROBLEM_H

#include <vector>
#include <string>
#include <filesystem>

using std::vector;
using std::string;
using std::filesystem::path;

enum class ProblemForm {
    NONE,
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
    static string form_to_string(ProblemForm form);

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


    auto getM1() const {
        vector<std::size_t> res;
        for(int i = 0; i < constraint_type_.size(); ++i) {
            if(constraint_type_[i] == ">=") {
                res.push_back(i);
            }
        }
        return res;
    }

    auto getM2() const {
        vector<std::size_t> res;
        for(int i = 0; i < constraint_type_.size(); ++i) {
            if(constraint_type_[i] == "=") {
                res.push_back(i);
            }
        }
        return res;
    }

    auto getN1() const { 
        vector<std::size_t> res;
        for(int i = 0; i < var_constraints_.size(); ++i) {
            if(var_constraints_[i] == ">=0") {
                res.push_back(i);
            }
        }
        return res;
    }

    auto getN2() const {
        vector<std::size_t> res;
        for(int i = 0; i < var_constraints_.size(); ++i) {
            if(var_constraints_[i] != ">=0") {
                res.push_back(i);
            }
        }
        return res;
    }

    ProblemForm getForm() const;

    void print(const string& title) const;
};



#endif
