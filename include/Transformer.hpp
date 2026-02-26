#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include "LinearProblem.hpp"

class Transformer {
private:
    static LinearProblem to_general_form(const LinearProblem& lp);
    static LinearProblem to_canonical_form(const LinearProblem& lp);
    static LinearProblem to_symmetrical_form(const LinearProblem& lp);
public:
    static LinearProblem to_general(const LinearProblem& lp);
    static LinearProblem to_canonical(const LinearProblem& lp);
    static LinearProblem to_symmetrical(const LinearProblem& lp);
};



#endif
