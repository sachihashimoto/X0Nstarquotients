 // load "QuadraticPoints/models_and_maps.m";
 load "J0wplusminus.m";

function XZeroNstar_nonhyperelliptic(N)
  number_of_terms := 20;
  //S := CuspForms(N);
  g := GenusX0NStar(N); //Dimension(S);
  //printf "g = %o\n", g;
  S, ALs := all_diag_basis(N);
  Sstar := [S[i] : i in [1..#S] | forall{w : w in ALs | w[i,i] eq +1}];
  assert #Sstar eq g;
  //bas2 := qExpansionBasis(CuspForms(N), number_of_terms);
  //print Type(bas), Type(bas[1]), Type(bas2), Type(bas2[1]);
  repeat
    number_of_terms +:= 10;
    //printf "Computing q-expansion basis of cusp forms of level %o of dimension %o with %o terms.\n", N, g, number_of_terms;
    //bas := qExpansionBasis(S, number_of_terms);
    bas := [qExpansion(f, number_of_terms) : f in Sstar];
    // determine equations for canonical model in P^{g-1}
    Pg1<[z]> := PolynomialRing(Rationals(), g);
    // intersection of hypersurfaces of degree d
    if g eq 3 then
      d:=4;
    elif g eq 4 or g eq 5 then
      d:=3;
    else
      d := 2;
    end if;
    //printf "Finding equation of model as intersection of %o.\n", d;
    monsd := MonomialsOfDegree(Pg1, d);
    mat := Matrix([[Coefficient(m, j) : j in [2..number_of_terms]] where m := Evaluate(mm, bas) : mm in monsd]);
    kermat := KernelMatrix(mat);
    //printf "Kernel matrix has %o rows and %o columns.\n", Nrows(kermat), Ncols(kermat);
    equations := [&+[kermat[i,j]*monsd[j] : j in [1..#monsd]] : i in [1..Nrows(kermat)]];
    Pws<q> := Universe(bas);
    X0_N_Scheme := Scheme(ProjectiveSpace(Pg1), equations);
    if Dimension(X0_N_Scheme) ne 1 then
      continue;
    end if;
    X0_N := Curve(X0_N_Scheme);
    // TODO: If this fails to often, X_0(N) is probably not canonically embedded.
  until assigned X0_N and Genus(X0_N) eq g;
  return X0_N;
end function;

function IsHyperellipticX0Nstar(N)
  if N in [1,2,3,4,5,6,7,8,9,10,12,13,16,18,25] cat  [11,14,15,17,19,20,21,24,27,32,36,49] cat [22,23,26,28,29,30,31,33,35,37,39,40,41,46,47,48,50,59,71] then
        return true;
    end if;
    if GenusX0NStar(N) le 2 or N in [136, 171, 176, 207, 252, 279, 315] then 
        return true;
    end if;
  return false;
end function;

function ratfun(target, basis, degree, R)
  // target = q-series of target function (e.g., j)
  // basis = [b1,...,bn] basis of space of cusp forms
  // degree = degree of numerator and denominator polynomial

  // try to write target as poly(basis)/poly(basis) with two
  // homogeneous polynomials of degree degree
  //P := PolynomialRing(Rationals(), #basis);
  mons := MonomialsOfDegree(R, degree);
  printf "Evaluating %o monomials of degree %o on basis of cardinality %o.\n", #mons, degree, #basis;
  evmons := [Evaluate(m, basis) : m in mons];
  evmonst := [-target*em : em in evmons];
  min := Min([Valuation(e) : e in evmons cat evmonst]);
  max := Min([AbsolutePrecision(e) : e in evmons cat evmonst]) - 1;
  mat := Matrix([[Coefficient(e, j) : j in [min..max]]
                   : e in evmons cat evmonst]);
  printf "Computing kernel matrix of matrix with %o rows and %o columns.\n", Nrows(mat), Ncols(mat);
  kermat := KernelMatrix(mat);
  printf "Kernel matrix with %o rows and %o columns.\n", Nrows(kermat), Ncols(kermat);
  result := [<&+[kermat[i,j]*mons[j] : j in [1..#mons]],
              &+[kermat[i,#mons+j]*mons[j] : j in [1..#mons]]> : i in [1..Nrows(kermat)]];
  result := [pair : pair in result | Valuation(Evaluate(pair[1], basis)) lt max];
  return result;
end function;

//function XZeroNstar_hyperelliptic(N)
/*  number_of_terms := 100;
 
 /*
  M := CuspForms(N);
  wds := [AtkinLehnerOperator(M, d) : d in Divisors(N) | d gt 1 and GCD(d, ExactQuotient(N,d)) eq 1];
  diag, D := Diagonalization(wds);
  i_s := [i : i in [1..Nrows(diag[1])] | forall{wd : wd in diag | wd[i,i] eq +1}];
  qbas := Basis(M, 2*number_of_terms);
  bas := [&+[D[i,j] * qbas[j] : j in [1..#qbas]] : i in i_s]; // TODO: take the right basis
  fs := [f : f in bas | Valuation(f) eq 1];
  f1 := fs[1];
  f2 := fs[2];
  Qq<q> := Parent(f1);
  g := #bas;
  // compute x, y in terms of f_1, f_2
  xq := f2/f1;
  //xq := Universe(bas)!xq;
  xq := Qq!xq;
  //yq := Parent(xq)!(Parent(xq).1*Derivative(xq)/bas[1]);
  yq := Qq!(q*Derivative(xq)/f1);*/
  // find a hyperelliptic equation of degree 2g + 2
 /* mat := Matrix([[Coefficient(xn, j) : j in [0..number_of_terms]] where xn := xq^n : n in [0..2*g+2]]);
  P<x> := PolynomialRing(Rationals());
  f1 := P!Eltseq(Solution(mat, Vector([Coefficient(y2, j) : j in [0..number_of_terms]]) where y2 := yq^2));
  C1 := HyperellipticCurve(f1);
  assert Conductor(C1) eq N^g;
  assert Genus(C1) eq g;
  // TODO: repeat with higher number_of_terms if output not correct
  return C1;*/
/*
  // compute the j-invariant morphism j_N: X -> P^1 as a function in xq, yq
  Pws<q> := Universe(bas);
  jq := jInvariant(q + O(q^(2*number_of_terms)));
  // f_1(q) &= q \cdot \frac{x'(q)}{y(q)}
  // f_2(q) &= x(q) \cdot f_1(q).

  degree := 0;
  repeat
    degree +:= 1;
    time jN := ratfun(jq, bas, degree);
  until not IsEmpty(jN);
  num, den := Explode(jN[1]);
  //num := Pg1!num; den := Pg1!den;
  return C1, num, den;//map<X0_N -> ProjectiveSpace(Rationals(),1) | [num, den]>;*/
//end function; 

function XZeroNstar(N)
    if IsHyperellipticX0Nstar(N) then 
      wds := [[d : d in Divisors(N) | Gcd(d, ExactQuotient(N,d)) eq 1]];
      X, ws, pairs, NB, cusp := eqs_quos(N, wds);
      X := pairs[1][1];
      return X;
      // return XZeroNstar_hyperelliptic(N);
    else
      // now X_0(N)^* has g > 2 and is not hyperelliptic
      X := XZeroNstar_nonhyperelliptic(N);
      return X;
    end if;
end function;
