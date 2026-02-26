#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "LinearProblem.hpp"

LinearProblem::LinearProblem(bool problem_type
                            ,vector<double> target
                            ,vector<vector<double>> constraints
                            ,vector<string> constraint_type
                            ,vector<double> rhs
                            ,vector<string> var_constraints)
        : problem_type_(problem_type)
        , target_(std::move(target))
        , constraints_(std::move(constraints))
        , constraint_type_(std::move(constraint_type))
        , rhs_(std::move(rhs))
        , var_constraints_(std::move(var_constraints)) {}

LinearProblem LinearProblem::read_from_file(const path& file) {
    if(!std::filesystem::exists(file)) {
        throw std::invalid_argument("File is not exists");
    }

    if(!std::filesystem::is_regular_file(file)) {
        throw std::invalid_argument("this file type is not supported");
    }

    std::ifstream stream(file);

    if(!stream.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    string line;
    int line_number = 1;

    if(!std::getline(stream, line) || (line != "min" && line != "max")) {
        throw std::runtime_error("Line 1: expected \"min\" or \"max\"");
    }
    line_number++;

    bool target_type = (line == "min"); 

    if(!std::getline(stream, line)) {
        throw std::runtime_error("Line 2: expected number of variables");
    }
    line_number++;

    int num;
    try {
        num = std::stoi(line);
        if(num <= 0) {
            throw std::invalid_argument("number of variables must be positive");
        }
    } catch(const std::exception& e) {
        throw std::runtime_error(string("Line 2: ") + e.what());
    }

    if(!std::getline(stream, line)) {
        throw std::runtime_error("Line 3: expected coefficients of target function");
    }
    line_number++;

    vector<double> target(num);
    std::istringstream iss_target(line);
    int counter = 0;

    while(counter < num && iss_target >> target[counter]) {
        counter++;
    }

    if(counter != num || !iss_target.eof()) {
        throw std::runtime_error("Line 3: wrong number of coefficients");
    }

    string extra;

    if(!std::getline(stream, line)) {
        throw std::runtime_error("Line 4: expected number of constraints");
    }
    line_number++;

    int constr_num;
    try {
        constr_num = std::stoi(line);
        if(constr_num <= 0) {
            throw std::invalid_argument("number of constraints must be positive");
        }
    } catch(const std::exception& e) {
        throw std::runtime_error(string("Line 4: ") + e.what());
    }

    vector<vector<double>> constraints(constr_num, vector<double>(num));
    vector<string> constraint_type(constr_num);
    vector<double> rhs(constr_num);

    for(int i = 0; i < constr_num; ++i) {
        if(!std::getline(stream, line)) {
            throw std::runtime_error("Line " + std::to_string(line_number) + ": expected constraint data");
        }

        std::istringstream iss_constraint(line);
        counter = 0;

        while(counter < num && iss_constraint >> constraints[i][counter]) {
            counter++;
        }

        if(counter != num) {
            throw std::runtime_error("Line " + std::to_string(line_number) +": wrong number of coefficients, expected ");
        }

        string type;
        if(!(iss_constraint >> type) ||
           (type != "=" && type != "<=" && type != ">=")) {
            throw std::runtime_error("Line " + std::to_string(line_number) + ": expected constraint type (=, <=, or >=)");
        }
        constraint_type[i] = type;

        if(!(iss_constraint >> rhs[i])) {
            throw std::runtime_error("Line " + std::to_string(line_number) + ": expected right-hand side value");
        }

        if(iss_constraint >> extra) {
            throw std::runtime_error("Line " + std::to_string(line_number) + ": extra data after right-hand side");
        }

        line_number++;
    }

    if(!std::getline(stream, line)) {
        throw std::runtime_error("Line " + std::to_string(line_number) + ": expected variable constraints");
    }

    vector<string> var_constraint(num);  
    std::istringstream iss_var(line);
    counter = 0;

    while(counter < num && iss_var >> var_constraint[counter]) {
        if(var_constraint[counter] != "<=0"
           && var_constraint[counter] != ">=0"
           && var_constraint[counter] != "free") {
            throw std::runtime_error("Line " + std::to_string(line_number) + ": invalid variable constraint '");
        }
        counter++;
    }

    if(counter != num || !iss_var.eof()) {
        throw std::runtime_error("Line " + std::to_string(line_number) +
                                ": wrong number of variable constraints, expected ");
    }

    return LinearProblem(target_type, target, constraints,
                        constraint_type, rhs, var_constraint);
}

LinearProblem LinearProblem::read_form_console() {
    std::cout << "=== Ввод задачи линейного программирования ===" << std::endl;
        
        std::string line;
        bool target_type;
        
        std::cout << "Введите тип задачи (min/max): ";
        std::getline(std::cin, line);
        
        if(line != "min" && line != "max") {
            throw std::invalid_argument("Invalid problem type. Expected 'min' or 'max'");
        }

        target_type = line == "min";
        
        int num_variables;
        std::cout << "Введите количество переменных: ";
        std::cin >> num_variables;
        
        if(std::cin.fail() || num_variables <= 0) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            throw std::invalid_argument("Invalid number of variables. Must be positive integer");
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        
        std::vector<double> target(num_variables);
        std::cout << "Введите " << num_variables << " коэффициентов целевой функции: ";
        
        std::getline(std::cin, line);
        std::istringstream iss_target(line);
        
        int count = 0;
        while(count < num_variables && iss_target >> target[count]) {
            count++;
        }
        
        if(count != num_variables || !iss_target.eof()) {
            throw std::invalid_argument("Expected " + std::to_string(num_variables) + 
                                       " coefficients for objective function, got " + 
                                       std::to_string(count));
        }
        
        
        int num_constraints;
        std::cout << "Введите количество ограничений: ";
        std::cin >> num_constraints;
        
        if(std::cin.fail() || num_constraints < 0) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            throw std::invalid_argument("Invalid number of constraints. Must be non-negative integer");
        }
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        
        std::vector<std::vector<double>> constraints(num_constraints, 
                                                     std::vector<double>(num_variables));
        std::vector<std::string> constraint_type(num_constraints);
        std::vector<double> rhs(num_constraints);
        
        std::cout << "\nВведите ограничения в формате: a1 a2 ... a" << num_variables 
                  << " {=, <=, >=} b" << std::endl;
        std::cout << "Каждое ограничение на отдельной строке:" << std::endl;
        
        for(int i = 0; i < num_constraints; ++i) {
            std::cout << "Ограничение " << (i + 1) << ": ";
            std::getline(std::cin, line);
            
            if(line.empty()) {
                throw std::invalid_argument("Empty constraint line");
            }
            
            std::istringstream iss_constraint(line);
            
            count = 0;
            while(count < num_variables && iss_constraint >> constraints[i][count]) {
                count++;
            }
            
            if(count != num_variables) {
                throw std::invalid_argument("Constraint " + std::to_string(i + 1) + 
                                           ": expected " + std::to_string(num_variables) + 
                                           " coefficients, got " + std::to_string(count));
            }
            
            std::string type;
            if(!(iss_constraint >> type)) {
                throw std::invalid_argument("Constraint " + std::to_string(i + 1) + 
                                           ": missing constraint type");
            }
            
            if(type != "=" && type != "<=" && type != ">=") {
                throw std::invalid_argument("Constraint " + std::to_string(i + 1) + 
                                           ": invalid constraint type '" + type + 
                                           "'. Expected '=', '<=', or '>='");
            }
            constraint_type[i] = type;
            
            double b;
            if(!(iss_constraint >> b)) {
                throw std::invalid_argument("Constraint " + std::to_string(i + 1) + 
                                           ": missing right-hand side value");
            }
            rhs[i] = b;
            
            if(!iss_constraint.eof()) {
                throw std::invalid_argument("Constraint " + std::to_string(i + 1) + 
                                           ": extra data after right-hand side");
            }
        }
        
        std::vector<std::string> var_constraints(num_variables);
        std::cout << "\nВведите ограничения на знак переменных (" << num_variables << " шт.)" << std::endl;
        std::cout << "Формат: для каждой переменной укажите: >=0, <=0, или free" << std::endl;
        std::cout << "Введите через пробел: ";
        
        std::getline(std::cin, line);
        std::istringstream iss_var(line);
        
        count = 0;
        while(count < num_variables && iss_var >> var_constraints[count]) {
            if(var_constraints[count] != ">=0" && 
               var_constraints[count] != "<=0" && 
               var_constraints[count] != "free") {
                throw std::invalid_argument("Variable " + std::to_string(count + 1) + 
                                           ": invalid constraint '" + var_constraints[count] + 
                                           "'. Expected '>=0', '<=0', or 'free'");
            }
            count++;
        }
        
        if(count != num_variables || !iss_var.eof()) {
            throw std::invalid_argument("Expected " + std::to_string(num_variables) + 
                                       " variable constraints, got " + std::to_string(count));
        }
        
        
        
        return LinearProblem(target_type, target, constraints, 
                            constraint_type, rhs, var_constraints);


}

ProblemForm LinearProblem::getForm() const {
    
    bool all_equalities = std::ranges::all_of(constraint_type_, [](const auto& str) { return str == "="; });
    
    bool all_nonnegative = std::ranges::all_of(var_constraints_, [](const auto& str) { return str == ">=0"; });
    
    if (all_equalities && all_nonnegative) {
        return ProblemForm::CANONIC;
    }
    
    bool all_geq = std::ranges::all_of(constraint_type_, [](const auto& str) { return str == ">="; });
    
    if (all_geq && all_nonnegative) {
        return ProblemForm::SYMMETRIC;
    }
    
    bool has_only_allowed_constraints = std::ranges::all_of(constraint_type_, [](const auto& str) { return str == ">=" || str == "="; }); 
    
    bool valid_sign_constraints = std::ranges::all_of(var_constraints_, [](const auto& str) { return str == ">=0" || str == "free"; });
    
    if (has_only_allowed_constraints && valid_sign_constraints) {
        return ProblemForm::GENERAL;
    }
    
    return ProblemForm::NONE;
}

string LinearProblem::form_to_string(ProblemForm f) {
    switch(f) {
        case ProblemForm::GENERAL:
            return "General";
        case ProblemForm::SYMMETRIC:
            return "Symmetric";
        case ProblemForm::CANONIC:
            return "Canonic";
        case ProblemForm::NONE:
            return "Nane";
    }
}


void LinearProblem::print(const std::string& title) const {
    std::cout << "\n" << std::string(50, '=') << std::endl;
    std::cout << title << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    std::cout << "Тип задачи: " << (!problem_type_ ? "MAXIMIZATION" : "MINIMIZATION") << std::endl;
    
    std::cout << "\nЦелевая функция:\n  ";
    if (!problem_type_) std::cout << "min  ";
    else std::cout << "max  ";
    
    std::cout << "Z = ";
    for (size_t i = 0; i < target_.size(); ++i) {
        if (i > 0) {
            if (target_[i] >= 0) std::cout << " + ";
            else std::cout << " - ";
        } else if (target_[i] < 0) {
            std::cout << "-";
        }
        
        double abs_val = std::abs(target_[i]);
        if (abs_val != 1.0) {
            std::cout << abs_val;
        }
        std::cout << "x" << (i + 1);
    }
    std::cout << std::endl;
    
    std::cout << "\nОграничения:\n";
    for (size_t i = 0; i < constraints_.size(); ++i) {
        std::cout << "  ";
        
        for (size_t j = 0; j < constraints_[i].size(); ++j) {
            if (j > 0) {
                if (constraints_[i][j] >= 0) std::cout << " + ";
                else std::cout << " - ";
            } else if (constraints_[i][j] < 0) {
                std::cout << "-";
            }
            
            double abs_val = std::abs(constraints_[i][j]);
            if (abs_val != 1.0) {
                std::cout << abs_val;
            }
            std::cout << "x" << (j + 1);
        }
        
        std::cout << " " << constraint_type_[i] << " " << rhs_[i] << std::endl;
    }
    
    std::cout << "\nОграничения на знак переменных:\n  ";
    for (size_t i = 0; i < var_constraints_.size(); ++i) {
        std::cout << "x" << (i + 1) << " " << var_constraints_[i];
        if (i < var_constraints_.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
    std::cout << std::string(50, '=') << "\n" << std::endl;
}





