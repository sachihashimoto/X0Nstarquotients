//276 star with MWS
N := 276;
SetLogFile("X0276star2.log");
HallDivs := [d : d in Divisors(N) | d ne 1 and GCD(d,N div d) eq 1];
X := SimplifiedModel(X0NQuotient(N,HallDivs));
J := Jacobian(X);
A, phi, boo1, boo2, rkbd := MordellWeilGroupGenus2(J);
n := #Generators(A);
assert Order(A.n) eq 0;
ans := Chabauty(phi(A.n));
printf "Number of rational points on X0(%o)star = %o\n", N, #ans;

N := 198;

SetLogFile("X0198star2.log");
HallDivs := [d : d in Divisors(N) | d ne 1 and GCD(d,N div d) eq 1];
X := SimplifiedModel(X0NQuotient(N,HallDivs));
J := Jacobian(X);
A, phi, boo1, boo2, rkbd := MordellWeilGroupGenus2(J);
n := #Generators(A);
assert Order(A.n) eq 0;
ans := Chabauty(phi(A.n));
printf "Number of rational points on X0(%o)star = %o\n", N, #ans;
