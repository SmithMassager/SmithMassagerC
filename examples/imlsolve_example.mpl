randomize():
with(LinearAlgebra):
with(ExternalCalling):

imlSolve := DefineExternal("imlSolve_maple", "../lib/libhnfproj.so"):

A := Matrix([[1,2], [0, 2]]);

b := Matrix([[1, 1],[2, 1]]);

imlSolve(A, b);
# Returns numerator matrix and a single integer reprsenting the denom
