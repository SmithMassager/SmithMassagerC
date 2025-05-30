macro(DefineExternal=ExternalCalling:-DefineExternal);

uniCertHelper := DefineExternal("unicert_maple", "lib/libhnfproj.so");
imlSolveHelper := DefineExternal("imlSolve_maple", "lib/libhnfproj.so");
highOrderResidue := DefineExternal("highOrderResidue_maple", "lib/libhnfproj.so");
imlMult := DefineExternal("imlMultiply_maple", "lib/libhnfproj.so");
add_maple := DefineExternal("add_maple", "lib/libhnfproj.so");
fastIntCertC := DefineExternal("fastIntCert_maple", "lib/libhnfproj.so");
iherm := DefineExternal("iherm_maple", "lib/libhnfproj.so");
# TODO: combine plMultiply and colplMultiply into one and similarly for the other ones.
plMultiply := DefineExternal("plMultiply_maple", "lib/libhnfproj.so");
colplMultiply := DefineExternal("colplMultiply_maple", "lib/libhnfproj.so");
cmodMulPL := DefineExternal("cmodMulviaPL_maple", "lib/libhnfproj.so");
cmodMulColPL := DefineExternal("cmodMulviaColPL_maple", "lib/libhnfproj.so");
cmodMulColPLMpz := DefineExternal("cmodMulviaPLmpz_maple", "lib/libhnfproj.so");
fmpzMult := DefineExternal("fmpzMult_maple", "lib/libhnfproj.so");
#mat_add := define_external('mat_add', LIB="lib/libhnfproj.so", 'A'::ARRAY(datatype=integer[8]), 'B'::ARRAY(datatype=integer[8]), 'C'::ARRAY(datatype=integer[8]), 'n'::integer[4], 'm'::integer[4]);
