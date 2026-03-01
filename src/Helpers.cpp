#include "Helpers.hpp"

void print_val_with_err(const double val1, const double val2) {
    double err = std::abs(val1 - val2) / 2;
    double avg = (val1 + val2) / 2;

    const double exponent = std::floor(std::log10(err));
    const auto precision = std::clamp(-exponent, 4.0, 10.0);
    const double factor = std::pow(10.0, -precision - 1); // Сохраняем 2 значащие цифры
    err = std::round(err / factor) * factor;
    avg = std::round(avg / factor) * factor;
    std::cout << std::fixed << std::setprecision(precision) << avg << " +- " << err << std::endl;
    const double rel_err = (err / std::abs(avg)) * 100.0;
    std::cout << "Real error: " << std::fixed << std::setprecision(2) << rel_err << " %" << std::endl;
}

std::vector<double> restore_original_solution(const LinearProblem& original_lp,
                                                     const std::vector<double>& canonical_solution) {
    const int n_orig = original_lp.num_variables();
    std::vector<double> original_solution(n_orig, 0.0);

    std::vector<int> free_neg_indices;
    for (int i = 0; i < n_orig; ++i) {
        if (original_lp.var_constaints()[i] == "free") {
            free_neg_indices.push_back(n_orig + static_cast<int>(free_neg_indices.size()));
        }
    }

    int free_counter = 0;
    for (int i = 0; i < n_orig; ++i) {
        if (original_lp.var_constaints()[i] == "free") {
            double pos_part = (i < static_cast<int>(canonical_solution.size())) ? canonical_solution[i] : 0.0;
            double neg_part = 0.0;

            if (free_counter < static_cast<int>(free_neg_indices.size()) &&
                free_neg_indices[free_counter] < static_cast<int>(canonical_solution.size())) {
                neg_part = canonical_solution[free_neg_indices[free_counter]];
            }

            original_solution[i] = pos_part - neg_part;
            free_counter++;
        } else {
            original_solution[i] = (i < static_cast<int>(canonical_solution.size())) ? canonical_solution[i] : 0.0;
        }
    }

    return original_solution;
}

