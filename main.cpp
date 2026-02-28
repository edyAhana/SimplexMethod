#include <iostream>

#include "LinearProblem.hpp"
#include "Transformer.hpp"
#include "Solver.hpp"

int main(int argc, char* argv[]) {
    LinearProblem lp = LinearProblem::read_from_file(argv[1]);
    auto gp =Transformer::to_general(lp); 
    auto cp =Transformer::to_canonical(lp); 
    auto sp =Transformer::to_symmetrical(lp); 

    lp.print("basic");
    Transformer::to_dual(lp).print("dual for basic");
    gp.print("general");
    Transformer::to_dual(gp).print("dual for general");
    sp.print("symmetric");
    Transformer::to_dual(sp).print("dual for symmetric");
    cp.print("canonic");
    Transformer::to_dual(cp).print("dual for canonic");

    Solver::solve(lp, true);
    Solver::solve(Transformer::to_dual(lp), true);

}
