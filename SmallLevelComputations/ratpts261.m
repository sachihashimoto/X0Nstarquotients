load "models_and_maps.m";

function MakeCharacter_87_a()
    N := 87;
    order := 1;
    char_gens := [59, 31];
    v := [1, 1];
    // chi(gens[i]) = zeta^v[i]
    assert UnitGenerators(DirichletGroup(N)) eq char_gens;
    F := CyclotomicField(order);
    chi := DirichletCharacterFromValuesOnUnitGenerators(DirichletGroup(N,F),[F|F.1^e:e in v]);
    return MinimalBaseRingCharacter(chi);
end function;

function all_diag_basis_AL(N, ds);
    C := CuspForms(N);
    g := Dimension(C);

    al_inds := [ m : m in Divisors(N) | GCD(m,N div m) eq 1 and m gt 1];
    al_invols := [AtkinLehnerOperator(C,d) : d in al_inds];

    T, new_als := simul_diag(al_invols);
    B := Basis(C);
    cleardenom := LCM([Denominator(x) : x in Eltseq(T)]);
    NB := [&+[cleardenom*T[i,j]*B[j] : j in [1..g]] : i in [1..g]];

    function is_invariant(i, d)
        return new_als[Index(al_inds, d)][i,i] eq +1;
    end function;
    invariant_indices := [i : i -> f in NB | forall{d : d in ds | is_invariant(i, d)}];
    assert #invariant_indices eq genus_quo(N, ds);
    NBwd := NB[invariant_indices];
    return NBwd;
end function;

function ConstructHyperellCurve(V :prec := 800)
    R<x> := PolynomialRing(Rationals());
    f1,f2 := Explode(qExpansionBasis(V, prec));
    x := f1/f2;
    K<q>:= Parent(x);
    y := Derivative(x)/f2*q;
    Monomials := [ 1, x , x^2 , x^3, x^4, x^5, x^6];
    Mat := [[Coefficient(Monomials[i], j): j in [0..400]] : i in [1..#Monomials]];
    Append(~Mat, [Coefficient(y^2, j) : j in [0..400]]);
    B := Kernel(Matrix(Mat));
    v := Basis(B)[1];
    R<t>:=PolynomialRing(Rationals());
    R!Eltseq([v[i]/v[8] : i in [1..7]]);
    f := -R!Eltseq([v[i]/v[8] : i in [1..7]]);
    H := HyperellipticCurve(f);
    return H, x, y;
end function;

function solve_x_in_canonical_ring(N, forms, xs)
    //solve for the sequence of things in x in the canonical ring as a rational function of forms
    d1 := 2;
    indexGam:=N*&*[Rationals() | 1+1/p : p in PrimeDivisors(N)];
    g := #forms;
    R<[a]> := PolynomialRing(Rationals(),g);

    done := false;
    while not done do
        printf "trying d1 = %o ...\n", d1;
        prec := Ceiling(indexGam* 2*d1/12);
        expanded_forms := [qExpansion(f, prec) : f in forms]; //all eigenforms
        monomials := MonomialsOfDegree(R,d1); //monomials of degree 2
        evals := [Evaluate(mon, expanded_forms) : mon in monomials]; //evaluate monomials of degree 2 in eigenforms
        evals_x := [];
        for x in xs do
            Append(~evals_x, [f*x : f in evals]); // multiply by x
        end for;
        allevals_x := &cat evals_x;
        prec := Minimum([Degree(f) : f in evals cat allevals_x ]);
        val := Minimum([Valuation(f) : f in evals cat allevals_x]);
        Pxs := [];
        Qxs := [];
        for i->x in xs do
            done := true;
            Mx := Matrix(Rationals(),[[Coefficient(f,i) : i in [val..prec]] : f in evals cat evals_x[i]]);
            if Rank(Mx) ge #monomials*2 then 
                printf "min(%o) >= %o\n", Rank(Mx), #monomials*2;
                d1 +:= 1;
                done := false;
                continue;
            end if;
            vx := ClearDenominator(Kernel(Mx).1)[1];
            Px := &+[monomials[i]*vx[i] : i in [1..#monomials]];
            Qx := &+[monomials[i]*vx[i+#monomials] : i in [1..#monomials]];
            Qx_eval := &+[Evaluate(monomials[i],expanded_forms)*vx[i+#monomials] : i in [1..#monomials]] ;

            if IsWeaklyZero(Qx_eval)  then
                printf "Some Qx = 0.\n";
                d1 +:= 1;
                done := false;
                continue;
            end if;
            Append(~Pxs, Px);
            Append(~Qxs, Qx);
        end for;

    end while;
    return Pxs, Qxs;

end function;

R<x> := PolynomialRing(Rationals());
chi := MakeCharacter_87_a();
Snew := NewSubspace(CuspidalSubspace(ModularSymbols(chi,2,0)));
V := Kernel([<2,R![-1, -1, 1]>],Snew);
H, x, y := ConstructHyperellCurve(V);

Hq := BaseChange(H, Parent(x));
Jq := Jacobian(Hq);
my_map := Hq![x,y];
M := 87;
N_over_M := 3;
d := 3;
plus_sign := false;
minus_eignenvals := [3];
deg_map := Hq![Evaluate(x, q^d), Evaluate(y,q^d)];
immersion := my_map - deg_map;
J := Jacobian(H);

N := 261;
ps := PrimeDivisors(N);
AL_N := [p^Valuation(N, p) : p in ps];
forms := all_diag_basis_AL(N, AL_N);

C := canonic(forms);


c0 := Coefficient(immersion[1],0);
c1 := Coefficient(immersion[1],1);
e0 := Coefficient(immersion[1],0);
e1 := Coefficient(immersion[1],1);

Ps, Qs := solve_x_in_canonical_ring(N, forms, [c0,c1, e0, e1]);

M<z> := PolynomialRing(Rationals());

P3<[x]> := Ambient(C);
P0x := Evaluate(Ps[1], [x[i]: i in [1..4]]);
Q0x := Evaluate(Qs[1], [x[i]: i in [1..4]]);
P1x := Evaluate(Ps[2], [x[i]: i in [1..4]]);
Q1x := Evaluate(Qs[2], [x[i]: i in [1..4]]);
R0x := Evaluate(Ps[3], [x[i]: i in [1..4]]);
R1x := Evaluate(Ps[4], [x[i]: i in [1..4]]);
S0x := Evaluate(Qs[3], [x[i]: i in [1..4]]);
S1x := Evaluate(Qs[4], [x[i]: i in [1..4]]);
FFC := FunctionField(C);

A1 := -FFC!(P1x/Q1x);
A0 := -FFC!(P0x/Q0x);
B1 := -FFC!(R1x/S1x);
B0 := -FFC!(R0x/S0x);

// M<z> := PolynomialRing(FFC);

// phi  := map<C -> J | [z^2 + A1*z +A0, B1*z + B0]>;

MW, mMW := MordellWeilGroupGenus2(J);
ptsMW := [p: p in MW];
pt := mMW(ptsMW[2]);
ptA1 := 1;
ptA0 := 1;
ptB1 := 1;
ptB0 := 2;
fib1 := [P1x + ptA1*Q1x, P0x +Q0x* ptA0, R1x +S1x*ptB1, R0x+S0x*ptB0];
alleqns := fib1 cat DefiningEquations(C);
RationalPoints(Scheme(P3, alleqns));

pt := mMW(ptsMW[5]);
ptA1 := 1;
ptA0 := 1;
ptB1 := -1;
ptB0 := -2;
fib1 := [P1x + ptA1*Q1x, P0x +Q0x* ptA0, R1x +S1x*ptB1, R0x+S0x*ptB0];
alleqns := fib1 cat DefiningEquations(C);
RationalPoints(Scheme(P3, alleqns));

//also need to consider the poles of the coordinates
poles := Q1x*Q0x*S1x*S0x;
alleqns := [poles] cat DefiningEquations(C);
RationalPoints(Scheme(P3,alleqns));


//but the other point on the MW group of J is not in the image of the map
