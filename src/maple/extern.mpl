macro(DefineExternal=ExternalCalling:-DefineExternal);

uniCertHelper := DefineExternal("unicert_maple", "lib/libhnfproj.so");
imlSolveHelper := DefineExternal("imlSolve_maple", "lib/libhnfproj.so");
highOrderResidue := DefineExternal("highOrderResidue_maple", "lib/libhnfproj.so");
imlMult := DefineExternal("imlMultiply_maple", "lib/libhnfproj.so");
add_maple := DefineExternal("add_maple", "lib/libhnfproj.so");
fastIntCertC := DefineExternal("fastIntCert_maple", "lib/libhnfproj.so");
iherm := DefineExternal("iherm_maple", "lib/libhnfproj.so");

# cmodMulPL, cmodMulColPL
## A, args[1] reprsenting a matrix with appropiate size.
## B, args[2] reprsenting a vector/matrix with appropiate size.
## F, args[3] optinal integer.
## Returns A.B cmod F
cmodMulPL := DefineExternal("cmodMulviaPL_maple", "lib/libhnfproj.so");
cmodMulColPL := DefineExternal("cmodMulviaColPL_maple", "lib/libhnfproj.so");
# fmpzMult
## A, args[1]
## B, args[2]
## Return A.B
fmpzMult := DefineExternal("fmpzMult_maple", "lib/libhnfproj.so");
