macro(DefineExternal=ExternalCalling:-DefineExternal);

imlSolveHelper := DefineExternal("imlSolve_maple", "../lib/libhnfproj.so"):

with(LinearAlgebra):

n := 3:

A := RandomMatrix(n,n):

b := RandomMatrix(n,10):

imlSolveHelper(A,b);

matMult := DefineExternal("matMult_maple", "../lib/libhnfproj.so"):
imlMult := DefineExternal("imlMultiply_maple", "../lib/libhnfproj.so"):

A[1,1] := 8888888888888888888888888888888888888;

mods(matMult(A, b, 8, 100), 8);

mods(100 * A . b, 30);
imlMult(A, b, 100, 30);
