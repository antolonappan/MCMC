# Project-code
Parameter estimation of cosmological model


This is a C++ Code written for fitting data using Markov Chain Monte Carlo methods and Bayesian statistics.In particular this is a multi threaded program, here two thread will generate two parameter chains concurrently and Gelman-Rubin diagnostics will be done in the third thread.

You can compile the code using command:g++ -pthread -std=c++11 <code file>

Use x.txt, y.txt, dy.txt with parabola.cpp

dataanalysis1.py & dataanalysis2.py is for data generated by project.cpp.Data for this code can be Hubble data.

