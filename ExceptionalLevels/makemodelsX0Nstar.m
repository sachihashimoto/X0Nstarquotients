load "X0Nstar-non-squarefree/magma/code/QuadraticPoints/models_and_maps.m";
load "X0Nstar-non-squarefree/magma/code/J0wplusminus.m";
load "X0Nstar-non-squarefree/magma/code/modelsX0Nstar.m";

SetLogFile("X0Nstar-models_chatelet.log");

for N in [40, 48, 72, 80, 88, 96, 100, 108, 112, 120, 135, 144, 147, 162, 196, 240] do
    if IsHyperellipticX0Nstar(N) then 
      wds := [[d : d in Divisors(N) | Gcd(d, ExactQuotient(N,d)) eq 1]];
      X, ws, pairs, NB, cusp := eqs_quos(N, wds);
      X := pairs[1][1];
      printf "N = %o: %o\n", N, X;
    else
      // now X_0(N)^* has g > 2 and is not hyperelliptic
      X := XZeroNstar(N);
      printf "N = %o: %o\n", N, X;
    end if;
    if N le 200 then 
      X2 := X0NQuotient(N, [d : d in Divisors(N) | Gcd(d, ExactQuotient(N,d)) eq 1 and d gt 1]);
      if not IsIsomorphic(X, X2) then 
        printf "WARNING: X0NQuotient is wrong!\n";
      end if;
    end if;
    g := GenusX0NStar(N);
    if g eq 0 then 
      printf "genus 0\n\n";
    elif g eq 1 then 
      print MordellWeilGroup(X), "\n\n";
    else 
      if IsHyperellipticX0Nstar(N) then
        pts := Points(X : Bound := 200);
      else
        X := Scheme(AmbientSpace(X), DefiningEquations(X));
        pts := PointSearch(X, 200);
      end if;
      printf "%o points: %o\n\n", #pts, pts;
    end if;
end for;
