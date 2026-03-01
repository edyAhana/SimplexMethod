cmake_minimum_required(VERSION 3.21)
project(lp)

set(CMAKE_CXX_STANDARD 23)

add_library(lp
    EigenHandler.cpp  
    EnumSolver.cpp    
    Helpers.cpp       
    LinearProblem.cpp 
    Solver.cpp        
    Transformer.cpp
)

target_include_directories(lp PUBLIC ../include /opt/homebrew/Cellar/eigen/5.0.1/include/eigen3/)

