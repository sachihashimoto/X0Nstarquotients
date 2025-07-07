load "../J0wplusminus.m";
SetDebugOnError(true);
SetLogFile("minimal.log");

function IsAlmostSquarefree(N)
    exponents := Sort([p[2] : p in Factorization(N)]);
    return #exponents ge 1 and #Indices(exponents, 1) eq (#exponents - 1) and exponents[#exponents] in {2,3};
end function;

function IsAlmostSquarefree2(N)
    exponents := Sort([p[2] : p in Factorization(N)]);
    return #exponents ge 1 and #Indices(exponents, 1) eq (#exponents - 1);
end function;

function IsExceptionalTuple(tup)
    if #tup eq 2 then
        if tup[1] eq 2 then
            return tup[2] in {3, 5, 7, 11, 23};
        elif tup[1] eq 3 then
            return tup[2] in {2, 5, 11};
        elif tup[1] eq 5 then
            return tup[2] in {2};
        elif tup[1] eq 7 then
            return tup[2] in {3};
        else
            return false;
        end if;
    elif #tup eq 3 then
        return tup[1] eq 2 and Set(tup[2..3]) eq {3,5};
    else
        return false;
    end if;
end function;

function IsExceptional(N : verbose := false)
    if IsSquarefree(N) then 
        return false;
    end if;
    if #PrimeDivisors(N) le 1 then // prime power?
        return false;
    end if;
    primesquaredivisors := [p[1] : p in Factorization(N) | p[2] gt 1];
    for p in primesquaredivisors do 
        if p notin {2,3,5,7,13} then 
            return false;
        end if;
        for M in [M : M in Divisors(N) | IsAlmostSquarefree(M) and IsDivisibleBy(M, p^2) and IsSquare(ExactQuotient(N, M))] do
            if verbose then 
                print Factorization(M);
            end if;
            for q in [l : l in PrimeDivisors(M) | l ne p] do
                if verbose then 
                    print [p,q];
                end if;
                if not IsExceptionalTuple([p,q]) then
                    return false;
                end if;
            end for;
            for q,r in [l : l in PrimeDivisors(M)] do
                if q eq r or q eq p or r eq p then 
                    continue;
                end if;
                if verbose then 
                    print [p,q,r];
                end if;
                if not IsExceptionalTuple([p,q,r]) then
                //    print [p,q,r], M;
                    return false;
                end if;
            end for;
        end for;
    end for;
    return true;
end function;


function IsExceptionalType1(N)
    // N = p^l * Q, l > 1, (p, Q) exceptional
    Nsqf, Nsq := SquarefreeFactorization(N);
    if Nsq eq 1 then 
        return false;
    end if;
    flag, p, l := IsPrimePower(Nsq);
    if not flag then 
        return false;
    end if;
    Q := ExactQuotient(Nsqf, p^Valuation(Nsqf, p));
    return IsExceptionalTuple([p] cat PrimeDivisors(Q));
end function;

function IsExceptionalType2(N)
    // N = p^m q^n, p,q small, m,n > 1. If m odd, (q,p) exceptional. If n odd, (p,q) exceptional.
    pqs := Factorization(N);
    if #pqs ne 2 then 
        return false;
    end if;
    p := pqs[1][1];
    m := pqs[1][2];
    q := pqs[2][1];
    n := pqs[2][2];
    if m le 1 or n le 1 then 
        return false;
    end if;
    if p notin {2,3,5,7,13} or q notin {2,3,5,7,13} then 
        return false;
    end if;
    if IsOdd(m) then 
        if not IsExceptionalTuple([q,p]) then 
            return false;
        end if;
    end if;
    if IsOdd(n) then 
        if not IsExceptionalTuple([p,q]) then 
            return false;
        end if;
    end if;
    return true;
end function;

function IsExceptionalType3(N)
    // N = p^m q^n r, m,n > 1 and (p,q,r) in {(2, 3, 5), (2, 3, 11), (2, 7, 3), (3, 5, 2)}. m even. n even in last three cases.
    pqs := Factorization(N);
    if #pqs ne 3 then 
        return false;
    end if;
    p := pqs[1][1];
    m := pqs[1][2];
    q := pqs[2][1];
    n := pqs[2][2];
    r := pqs[3][1];
    s := pqs[3][2];
    if p ne 2 then
        return false;
    end if;
    if [q,r] eq [3,5] then //in {[3,5],[3,11]} then 
        return (m gt 1 and n gt 1 and s eq 1 and IsEven(m)) or
               (m eq 1 and n gt 1 and s gt 1 and IsEven(n) and IsEven(s));
    end if;
    if [q,r] eq [3,11] then 
        return m gt 1 and n gt 1 and s eq 1 and IsEven(m) and IsEven(n);
    end if;
    if [r,q] eq [7,3] then 
        return m gt 1 and s gt 1 and n eq 1 and IsEven(m) and IsEven(s);
    end if;
    return false;
end function;

function IsExceptionalType4(N)
    // N = R^2 with PrimeDivisors(R) small and # >= 2
    flag, R := IsSquare(N);
    if not flag then 
        return false;
    end if;
    primedivisors := Set(PrimeDivisors(R));
    return #primedivisors gt 1 and primedivisors subset {2,3,5,7,13};
end function;

function IsExceptionalType5(N)
    // N = 2^m 3^n 5^s, m odd, n, s even.
    if PrimeDivisors(N) eq [2,3,5] then 
        pqr := Factorization(N);
        return (IsOdd(pqr[1][2]) and pqr[1][2] ge 3 and IsEven(pqr[2][2]) and IsEven(pqr[3][2])) or
           (IsEven(pqr[1][2]) and IsEven(pqr[2][2])  and IsOdd(pqr[3][2]) and pqr[3][2] ge 3);
    elif PrimeDivisors(N) eq [2,3,7] then
        pqr := Factorization(N);
        return (IsEven(pqr[1][2]) and IsOdd(pqr[2][2]) and pqr[2][2] ge 3 and IsEven(pqr[3][2]));
    end if;
    return false;
end function;

// tests all the cases
// should be equivalent to IsExceptional
function IsExceptional2(N)
    return IsExceptionalType1(N) or IsExceptionalType2(N) or IsExceptionalType3(N) or IsExceptionalType4(N) or IsExceptionalType5(N);
end function;


for N in [2..10000] do
    if IsExceptional(N) eq IsExceptional2(N) then
        continue;
    end if;
    printf "N = %o = %o: %o ? %o = %o, %o, %o, %o, %o\n", N, Factorization(N), IsExceptional(N), IsExceptional2(N),
        IsExceptionalType1(N), IsExceptionalType2(N), IsExceptionalType3(N), IsExceptionalType4(N), IsExceptionalType5(N);
end for;



////////////////////////////////////
/// find minimal exceptional levels


//g0Nstarnull := [N : N in [2..119] | GenusX0NStar(N) eq 0]; // we know 199 is an upper bound
g0Nstarnull := [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 
39, 41, 42, 44, 45, 46, 47, 49, 50, 51, 54, 55, 56, 59, 60, 62, 66, 69, 70, 71, 
78, 87, 92, 94, 95, 105, 110, 119 ];
//printf "g_0(N)^* = 0 for N in %o\n", g0Nstarnull;

minimal_levels := {};


// exceptional levels of type 1
for p in [2,3,5,7,13] do
    if p eq 2 then 
        qs := [3, 5, 7, 11, 23];
    elif p eq 3 then
        qs := [2, 5, 11];
    elif p eq 5 then
        qs := [2];
    elif p eq 7 then
        qs := [3];
    else
        break;
    end if;
    for q in qs do
        l := 2;
        while true do
            N := p^l * q;
            assert IsExceptional(N);
            if N notin g0Nstarnull then
                //print N, Factorization(N);
                Include(~minimal_levels, N);
                break;
            end if;
            l +:= 2;
        end while;

        l := 3;
        while true do
            N := p^l * q;
            assert IsExceptional(N);
            if N notin g0Nstarnull then
                //print N, Factorization(N);
                Include(~minimal_levels, N);
                break;
            end if;
            l +:= 2;
        end while;
    end for;
end for;
printf "minimal levels 1.1, %o\n", minimal_levels;

// case 1.2 (exceptional triples, Q = 15)
p := 2;
Q := 15;
l := 2;
while true do
    N := p^l * Q;
    assert IsExceptional(N);
    if N notin g0Nstarnull then
        //print N, Factorization(N);
        Include(~minimal_levels, N);
        break;
    end if;
    l +:= 2;
end while;
printf "minimal levels 1.2, %o\n", minimal_levels;
l := 3;
while true do
    N := p^l * Q;
    assert IsExceptional(N);
    if N notin g0Nstarnull then
        //print N, Factorization(N);
        Include(~minimal_levels, N);
        break;
    end if;
    l +:= 2;
end while;
printf "minimal levels 1.2, %o\n", minimal_levels;

// type 2: N = p^m q^n, p,q small, m,n > 1. If m odd, (q,p) exceptional. If n odd, (p,q) exceptional.
for p in [2,3,5,7,13] do
    for q in [2,3,5,7,13] do
        if q eq p then 
            continue;
        end if;
        m := 2;
        n := 2;
        while true do
            N := p^m * q^n;
            //print N, Factorization(N);
            assert IsExceptional(N);
            if N notin g0Nstarnull then
                Include(~minimal_levels, N);
                break;
            end if;
            m +:= 2;
        end while;
    end for;
end for;
printf "minimal levels 2, %o\n", minimal_levels;
for p in [2,3,5,7,13] do
    for q in [2,3,5,7,13] do
        if q eq p then 
            continue;
        end if;
        m := 3;
        n := 2;
        if not IsExceptionalTuple([q,p]) then 
            continue;
        end if;
        while true do
            N := p^m * q^n;
            //print N, Factorization(N);
            assert IsExceptional(N);
            if N notin g0Nstarnull then
                Include(~minimal_levels, N);
                break;
            end if;
            m +:= 2;
        end while;
    end for;
end for;
printf "minimal levels 2, %o\n", minimal_levels;
for p in [2,3,5,7,13] do
    for q in [2,3,5,7,13] do
        if q eq p then 
            continue;
        end if;
        m := 2;
        n := 3;
        if not IsExceptionalTuple([p,q]) then 
            continue;
        end if;
        while true do
            N := p^m * q^n;
            //print N, Factorization(N);
            assert IsExceptional(N);
            if N notin g0Nstarnull then
                Include(~minimal_levels, N);
                break;
            end if;
            n +:= 2;
        end while;
    end for;
end for;
printf "minimal levels 2, %o\n", minimal_levels;
for p in [2,3,5,7,13] do
    for q in [2,3,5,7,13] do
        if q eq p then 
            continue;
        end if;
        m := 3;
        n := 3;
        if not IsExceptionalTuple([p,q]) or not IsExceptionalTuple([q,p]) then 
            continue;
        end if;
        while true do
            N := p^m * q^n;
            //print N, Factorization(N);
            assert IsExceptional(N);
            if N notin g0Nstarnull then
                Include(~minimal_levels, N);
                break;
            end if;
            n +:= 2;
        end while;
    end for;
end for;
printf "minimal levels 2, %o\n", minimal_levels;

// type 3: N = p^m q^n r, m,n > 1 and (p,q,r) in {(2, 3, 5), (2, 3, 11), (2, 7, 3), (3, 5, 2)}. m even, n even in last three cases.
for pqr in [[2,3,5], [2,3,11], [2,7,3], [3,5,2]] do 
    p,q,r := Explode(pqr);
    m := 2;
    n := 2;
    while true do
        N := p^m * q^n * r;
        //print N, Factorization(N);
        assert IsExceptional(N);
        if N notin g0Nstarnull then
            Include(~minimal_levels, N);
            break;
        end if;
        m +:= 2;
    end while;
end for;
printf "minimal levels 3, %o\n", minimal_levels;
for pqr in [[2,3,5], [2,3,11], [2,7,3], [3,5,2]] do 
    p,q,r := Explode(pqr);
    m := 2;
    n := 2;
    while true do
        N := p^m * q^n * r;
        //print N, Factorization(N);
        assert IsExceptional(N);
        if N notin g0Nstarnull then
            Include(~minimal_levels, N);
            break;
        end if;
        n +:= 2;
    end while;
end for;
printf "minimal levels 3, %o\n", minimal_levels;
// *not* in last three cases
for pqr in [[2,3,5]] do //, [2,3,11]] do 
    p,q,r := Explode(pqr);
    m := 2;
    n := 3;
    while true do
        N := p^m * q^n * r;
        //print N, Factorization(N);
        assert IsExceptional(N);
        if N notin g0Nstarnull then
            Include(~minimal_levels, N);
            break;
        end if;
        n +:= 2;
    end while;
end for;
printf "minimal levels 3, %o\n", minimal_levels;

// type 4: N = R^2 with PrimeDivisors(R) small.
for pps in Subsets({2,3,5,7,13}) do
    if #pps lt 2 then 
        continue;
    end if;
    ps := Setseq(pps);
    for i in [1..#ps] do
        ns := [2 : j in [1..#ps]];
        while true do
            N := &*[ps[j]^ns[j] : j in [1..#ps]];
            //print N, Factorization(N);
            assert IsExceptional(N);
            if N notin g0Nstarnull then
                Include(~minimal_levels, N);
                break;
            end if;
            ns[i] +:= 2; 
        end while;
    end for;
end for;

// type 5: N = 2^m 3^n 5^s, m odd, n, s even.
m := 1;
n := 2;
s := 2;
while true do
    N := 2^m * 3^n * 5^s;
    //print N, Factorization(N);
    assert IsExceptional(N);
    if N notin g0Nstarnull then
        Include(~minimal_levels, N);
        break;
    end if;
    m +:= 2;
end while;
printf "minimal levels 5, %o\n", minimal_levels;
m := 1;
n := 2;
s := 2;
while true do
    N := 2^m * 3^n * 5^s;
    //print N, Factorization(N);
    assert IsExceptional(N);
    if N notin g0Nstarnull then
        Include(~minimal_levels, N);
        break;
    end if;
    n +:= 2;
end while;
printf "minimal levels 5, %o\n", minimal_levels;
m := 1;
n := 2;
s := 2;
while true do
    N := 2^m * 3^n * 5^s;
    //print N, Factorization(N);
    assert IsExceptional(N);
    if N notin g0Nstarnull then
        Include(~minimal_levels, N);
        break;
    end if;
    s +:= 2;
end while;
printf "minimal levels 5, %o\n", minimal_levels;

function ExceptionalQuotients(N)
    quotients := [];
    for p in PrimeDivisors(N) do
        if not IsDivisibleBy(N, p^2) then
            continue;
        end if;
        M := ExactQuotient(N, p^2);
        if IsExceptional(M) then 
            if M notin g0Nstarnull then
                Append(~quotients, M);
            end if;
        end if;
    end for;
    further_quotients := &cat[ExceptionalQuotients(M) : M in quotients];
    if #further_quotients eq 0 then 
        return [N];
    end if;
    return further_quotients;
end function;

minimal_levels2 := {};
for N in minimal_levels do
    //if exists{p : p in PrimeDivisors(N) | IsDivisibleBy(N, p^2) and ExactQuotient(N, p^2) in g0Nstarnull} then
    //    Include(~minimal_levels2, N);
    //end if;
    Include(~minimal_levels2, Minimum(ExceptionalQuotients(N)));
end for;
minimal_levels := {};
for N in minimal_levels2 do
    if IsDivisibleBy(N, 13^2) and IsSquare(ExactQuotient(N, 13^2)) then
        Include(~minimal_levels, 13^2);
        continue;
    end if;
    Include(~minimal_levels, N);
end for;

// filter p^k * squarefree
minimal_levels2 := {};
for N in minimal_levels do
    if not IsAlmostSquarefree2(N) then
        Include(~minimal_levels2, N);
    end if;
end for;
print "minimal levels: ", minimal_levels,
        "minimal levels not almost squarefree: ", minimal_levels2,
        "minimal levels almost squarefree: ", minimal_levels diff minimal_levels2;
//minimal_levels := minimal_levels2;


minimal_levels := Sort(Setseq(minimal_levels));
//printf "reduced set of minimal exceptional levels:\n%o\n", minimal_levels;
printf "factorizations:\n%o\n", [Factorization(N) : N in minimal_levels];
//printf "<N, g_0(N)^*>:\n%o\n", [<N, GenusX0NStar(N)> : N in minimal_levels];

//printf "g_0(N)^*:\n";


function cusps(N)
    //cs := {<a,b> : a,b in [0..N] | (b eq 0 or (IsDivisibleBy(N, b) and a lt Gcd(b, ExactQuotient(N, b)))) and Gcd(a,b) eq 1};
    cs1 := {};
    for d in Divisors(N) do
        w := Gcd(d, ExactQuotient(N, d));
        if w eq 1 then 
            if d eq 1 then 
                Include(~cs1, <1,0>);
            elif d eq N then 
                Include(~cs1, <0,1>);
            else
                Include(~cs1, <1,d>);
            end if;
        else
            for a in [1..w] do 
                a1 := a;
                if Gcd(a1, w) eq 1 then
                    while Gcd(a1, ExactQuotient(d, w)) ne 1 do 
                        a1 +:= w;
                    end while;
                    Include(~cs1, <a1,d>);
                end if;
            end for;
        end if;
    end for;
    //print cs, cs1;
    assert #cs1 eq &+[EulerPhi(Gcd(d, ExactQuotient(N,d))) : d in Divisors(N)];
    return cs1;
end function;

function is_equivalent(c1, c2, N)
    if Gcd(c1[2], N) ne Gcd(c2[2], N) then 
        return false;
    end if;
    d := Gcd(c1[2], N);
    ZN := Integers(Gcd(d, ExactQuotient(N,d)));
    if ZN!(c1[1] * ExactQuotient(c1[2], d)) eq ZN!(c2[1] * ExactQuotient(c2[2], d)) then 
        return true;
    end if;
    return false;
end function;

function reduce_cusp(a, b, N, M)
    if <a,b> eq <1,0> then 
        return <a,b>;
    end if;
    cs := {};
    for c in cusps(M) do
        if is_equivalent(c, <a,b>, M) then
            Include(~cs, c);
        end if;
    end for;
    assert #cs eq 1;
    return Setseq(cs)[1];
end function;

function number_of_Galois_orbits_of_cusps_on_X0NStar(N)
    return &*[Integers()| Ceiling((pn[2] + 1)/2) : pn in Factorization(N)];
end function;

function is_Hall_divisor(N, d)
    return Gcd(d, ExactQuotient(N, d)) eq 1;
end function;

function HallDivisors(N) 
    return [q : q in Divisors(N) | is_Hall_divisor(N, q)];
end function;

function wQN(N, q)
    assert is_Hall_divisor(N, q);
    d, x, y:= XGCD(q, -ExactQuotient(N,q));
    assert d eq 1;
    g := Matrix(Integers(), 2,2, [q*x, y, N, q]);
    assert Determinant(g) eq q;
    return g;
end function;

function act_AtkinLehner_prime(c, p, N)
    assert IsPrime(p);
    alpha := Valuation(N, p);
    Na := ExactQuotient(N, p^alpha);
    wcs := {};
    a, b := Explode(c);
    for aprime,bprime in [0..N] do
        if b ne 0 then
            beta := Valuation(b, p);
            b0 := ExactQuotient(b, p^beta);
        else 
            beta := alpha;
            b0 := 0;
        end if;
        if Integers(Na)!aprime eq Integers(Na)!a and Integers(Na)!bprime eq Integers(Na)!b and
           Integers(p^alpha)!aprime eq Integers(p^alpha)!(-b0) and Integers(p^alpha)!bprime eq Integers(p^alpha)!(p^(alpha - beta) * a) then 
           Include(~wcs, <aprime,bprime>);
        end if;
    end for;
    wcs := {reduce_cusp(c[1], c[2], N, N) : c in wcs};
    assert #wcs eq 1;
    return Setseq(wcs)[1];
end function;

function act_AtkinLehner(c, d, N)
    assert is_Hall_divisor(N, d);
    wc := c;
    for p in PrimeDivisors(d) do 
        wc := act_AtkinLehner_prime(wc, p, N);
    end for;
    return wc;
end function;

function equivalence_classes_of_cusps_on_X0NStar(N)
    equivalence_classes := {{act_AtkinLehner(c, d, N) : d in HallDivisors(N)} : c in cusps(N)};
    assert &join equivalence_classes eq cusps(N);
    // are the equivalence classes pairwise disjoint?
    assert forall{class : class in equivalence_classes | forall{class2 : class2 in equivalence_classes diff {class} | class meet class2 eq {}}};
    return equivalence_classes;
end function;

// sanity check
for N in [2..100] do
    if IsSquarefree(N) then 
        assert #equivalence_classes_of_cusps_on_X0NStar(N) eq 1;
        /*if #equivalence_classes_of_cusps_on_X0NStar(N) ne 1 then 
            print N, number_of_Galois_orbits_of_cusps_on_X0NStar(N), #equivalence_classes_of_cusps_on_X0NStar(N), equivalence_classes_of_cusps_on_X0NStar(N);
        end if;*/
    end if;
    //assert number_of_Galois_orbits_of_cusps_on_X0NStar(N) eq #equivalence_classes_of_cusps_on_X0NStar(N);
    /*if number_of_Galois_orbits_of_cusps_on_X0NStar(N) ne #equivalence_classes_of_cusps_on_X0NStar(N) then 
        print N, number_of_Galois_orbits_of_cusps_on_X0NStar(N), #equivalence_classes_of_cusps_on_X0NStar(N), equivalence_classes_of_cusps_on_X0NStar(N);
    end if;*/
end for;

function act_character(c, d, N)
    return reduce_cusp(c[1] * d, c[2], N, N);
end function;

function act_Galois(c, N)
    return {act_character(c, d, N) : d in [d : d in [1..N-1] | Gcd(d,N) eq 1]};
end function;

function is_Galois_stable(c, N)
    return forall{d : d in c | act_Galois(d, N) subset c};
end function;

function Q_cusps_on_X0NStar(N)
    return {c : c in equivalence_classes_of_cusps_on_X0NStar(N) | is_Galois_stable(c, N)};
end function;

/*for N in [40, 48, 72, 80, 88, 96, 99, 100, 108, 112, 120, 135, 144, 147, 162, 176, 180, 184, 196, 200, 216, 224, 225, 240] do
    printf "N = %o: #C_0(N)^*(Q) = %o\n%o\n", N, #Q_cusps_on_X0NStar(N), Q_cusps_on_X0NStar(N);
end for;*/

function width(a, b, N);
  assert (b eq 0 and a eq 1) or (IsDivisibleBy(N, b) and Gcd(a,b) eq 1);
  return ExactQuotient(N, Gcd(N, b^2));
end function;

function is_unramified(a, b, N, M, d);
  assert (b eq 0 and a eq 1) or IsDivisibleBy(ExactQuotient(N, M), d);
  if <a,b> eq <1,0> then 
    return true;
  end if;
  assert IsDivisibleBy(N, b) and Gcd(a,b) eq 1;
  //if d eq 1 then 
  //  return ExactQuotient(width(a, b, N), width(reduce_cusp(a, b, N, M)[1], reduce_cusp(a, b, N, M)[2], M)) eq 1;
  //elif IsDivisibleBy(ExactQuotient(N,M), d) and width(a, b, N) eq 1 then
  //  return forall{p : p in PrimeDivisors(d) | Valuation(b,p) eq Valuation(N,p)/2 and Valuation(ExactQuotient(N,M), p) eq Valuation(d, p)};
  //else
  //  error "not implemented";
  //end if;
  // we use it only for width 1 cusps
  assert width(a, b, N) eq 1;
  return forall{p : p in PrimeDivisors(d) | Valuation(b,p) eq Valuation(N,p)/2 and Valuation(ExactQuotient(N,M), p) eq Valuation(d, p)};
end function;

function cusp_to_infty(a,b,N) // a matrix in Gamma_0(N) transporting a/b to 1/0
    //d, f, e := XGCD(a,-b*N);
    d, f, e := XGCD(a,-b);
    assert d eq 1;
    //g := GL(2,Integers())![f,-e,-b*N,a];
    g := GL(2,Integers())![f,-e,-b,a];
    assert Determinant(g) eq 1;
    //assert IsDivisibleBy(g[2,1], N);
    aprime := g[1,1] * a + g[1,2] * b;
    bprime := g[2,1] * a + g[2,2] * b;
    denom := Gcd(aprime, bprime);
    aprime := ExactQuotient(aprime, denom);
    bprime := ExactQuotient(bprime, denom);
    assert aprime eq 1 and bprime eq 0;
    return g;
end function;

function phiNM(M,q)
    wqM := &*[Integers()| p^Valuation(M,p) : p in PrimeDivisors(q)];
    assert is_Hall_divisor(M,wqM);
    //assert forall{p : p in PrimeDivisors(M) | IsDivisibleBy(wqM,p) eq IsDivisibleBy(M,p)};
    return wqM;
end function;

function root_of_unity(a,b,d,N,M,w_NM)
    assert is_Hall_divisor(ExactQuotient(N,M), d);
    assert is_unramified(a,b, N, M, d) and width(a,b,N) eq 1;
    dN := &*[Integers()| p^Valuation(N,p) : p in PrimeDivisors(d)];
    wdN := wQN(N, dN);
    aprimeprime := wdN[1,1] * a + wdN[1,2] * b;
    bprimeprime := wdN[2,1] * a + wdN[2,2] * b;
    //print a,b, wdN, aprime, bprime;
    denom := Gcd(aprimeprime, bprimeprime);
    aprimeprime := ExactQuotient(aprimeprime, denom);
    bprimeprime := ExactQuotient(bprimeprime, denom);
    //print aprime, bprime;
    aprime, bprime := Explode(reduce_cusp(aprimeprime, bprimeprime, N, N));
    //print aprime, bprime;
    assert bprime eq 0 or IsDivisibleBy(N, bprime);
  //  assert is_unramified(aprime,bprime, N, M, dM) and width(aprime,bprime,N) eq 1;
    //aprime, bprime := Explode(reduce_cusp(aprime, bprime, N, N));
    assert bprime eq b;
    assert is_unramified(aprime,b, N, M, d) and width(aprime,b,N) eq 1;
    gamma := cusp_to_infty(aprimeprime,bprimeprime,N);
    gammaprime := cusp_to_infty(aprime,b,N);
    // gammaprimeprime * a''/b''' = a'/b'
    gammaprimeprime := gammaprime^-1 * gamma;
    assert gammaprimeprime in Gamma0(N);
    assert gammaprimeprime[1,1] * aprimeprime + gammaprimeprime[1,2] * bprimeprime eq aprime and 
            gammaprimeprime[2,1] * aprimeprime + gammaprimeprime[2,2] * bprimeprime eq b;  
    wdN2 := gammaprimeprime * wdN;
    //assert wdN2[1,1] * a + wdN2[1,2] * b eq aprime and wdN2[2,1] * a + wdN2[2,2] * b eq b;     
    gamma := cusp_to_infty(a,b,N);
    gammaprime := cusp_to_infty(aprime,b,N);
    g := gammaprime * wdN2 * gamma^-1;
    //print g;
    assert g[2,1] eq 0 and g[1,1] eq g[2,2] and Gcd(g[1,1], g[1,2]) eq 1;
    assert Abs(g[1,1]) eq &*[Integers()| p^ExactQuotient(Valuation(N,p),2) : p in PrimeDivisors(d)];
    w := w_NM^ExactQuotient(ExactQuotient(N,M), Abs(g[1,1]));
    //K<w> := CyclotomicField(Abs(g[1,1]));
    assert w^Abs(g[1,1]) eq 1 and forall{e : e in Divisors(Abs(g[1,1])) | e eq Abs(g[1,1]) or w^e ne 1};
    return w^Abs(g[1,2]);
end function;

//for N in minimal_levels do
for N in [250,297,368,405,441,450,486,891,1029,1225,1250] do
    printf "N = %o = %o\n", N, Factorization(N);

    for M in Divisors(N) do 
        if M in {1,N} then 
            continue;
        end if;

        // hypothesis (*)
        //if not forall{p : p in {p : p in PrimeDivisors(N) | IsEven(Valuation(N,p))} | Valuation(M,p) le ExactQuotient(Valuation(N,p), 2)} then 
        if not forall{p : p in PrimeDivisors(N) | Valuation(M,p) le Ceiling(Valuation(N,p) / 2)} then 
            continue;
        end if;

        // hypothesis (RZQ)
        rank0new, rank0s, JMs := HasNewRank0Quotient(JZero(M));
        if rank0new then 
            function AtkinLehnerEigenvalue(JM, q)
                w := AtkinLehnerOperator(JM, q);
                charpol := CharacteristicPolynomial(w);
                roots := Roots(charpol);
                assert #roots eq 1;
                return roots[1][1];
            end function;
            //printf "with AL-eigenvalues %o.\n", [* [<q, AtkinLehnerEigenvalue(JM, q)> : q in HallDivisors(Level(JM))] : JM in JMs *];

            //printf "sums of roots of unity:\n";
            width1cusps := {c : c in cusps(N) | width(c[1], c[2], N) eq 1 and (c[2] ne 0 and IsDivisibleBy(N,c[2]))};

            for k -> JM in JMs do
            if not exists{q : q in HallDivisors(ExactQuotient(N, M)) | (q eq 1 or IsPrimePower(q)) and IsDivisibleBy(M, q) 
                    and AtkinLehnerEigenvalue(JM, phiNM(Level(JM), q)) eq -1} then 
                continue;
            end if;
            if IsDivisibleBy(N, 2) and IsEven(Valuation(N,2)) and AtkinLehnerEigenvalue(JM, 2^Valuation(Level(JM), 2)) ne -1 then
                continue;
            end if;
            printf "J_0(%o)^new has rank 0 quotients: %o\n", M, rank0s;
            //printf "rank 0 quotient #%o with AL-eigenvals %o:\n", k, [<q, AtkinLehnerEigenvalue(JM, q)> : q in HallDivisors(Level(JM))];
            ord_w := ExactQuotient(N, M);
            K<w> := CyclotomicField(ord_w);
            OK := Integers(K);
           // printf "w has degree %o\n", Degree(MinimalPolynomial(w));
            assert w^ord_w eq 1;
            lcm_n_c := 1;
            for c in width1cusps do 
                //if c ne <1,90> then continue; end if;
                Rc := {};
                for d in {d : d in Divisors(ExactQuotient(N, M))} do // TODO: "full valuations"
                  if is_Hall_divisor(ExactQuotient(N, M), d) and is_unramified(c[1], c[2], N, M, d) then
                    /*if not reduce_cusp(c[1], c[2], N, M) eq <1,0> then
                        print c, N, M, d;
                        assert false;
                    end if;*/
                    assert reduce_cusp(c[1], c[2], N, M) eq <1,0>;
                    //printf "unramified cusp %o for d = %o.\n", c, d;
                    Include(~Rc, d);
                  end if;
                end for;
                printf "R_c = %o\n", Rc;
                /*function possible_sums(summands)
                    if #summands eq 0 then
                        return {};
                    end if;
                    if #summands eq 1 then 
                        return summands[1];
                    else 
                        return {s1 + s2 : s1 in summands[1], s2 in possible_sums(summands[2..#summands])};
                    end if;
                end function;
                SignCharacters := [ func < n | &*[Integers()| ((-1)^(p in e select 1 else 0))^(IsDivisibleBy(n, p) select 1 else 0) : p in PrimeDivisors(ord_w)] > 
                    : e in Subsets(Set(PrimeDivisors(ord_w)))];
                for eps in SignCharacters do
                    if (eps(1) ne +1) or (eps(2) eq +1) then 
                        continue;
                    end if;
                    summands := [];
                    for d in Rc do
                        ord_zeta := &*[Integers()| p^ExactQuotient(Valuation(N, p), 2) : p in PrimeDivisors(d)];
                        zeta := w^(ExactQuotient(ord_w, ord_zeta));
                        assert zeta^ord_zeta eq 1;
                        print root_of_unity(c[1],c[2],d,N);
                        Append(~summands, {eps(d) * zeta^i : i in [1..ord_zeta] | Gcd(i, ord_zeta) eq 1});
                    end for;
                    printf "sums of roots of unity for c = %o:\n%o\n", c,
                     { Nm_s ne 0 select PrimeDivisors(Nm_s) else [0] where Nm_s := Integers()!Norm(s) : s in possible_sums(summands) };
                end for;*/
                //printf "signs = %o\n", [<d, phiNM(M,d), AtkinLehnerEigenvalue(JM, phiNM(M,d))> : d in Rc];
                roots_of_unity := [K| AtkinLehnerEigenvalue(JM, phiNM(M,d)) * root_of_unity(c[1],c[2],d,N,M, w) : d in Rc];
                //print [<e, Min([i : i in [1..10*ord_w] | e^i eq 1])> : e in roots_of_unity];
                S_c := &+roots_of_unity;
                //OK := Integers(K);
                Nm_s := Integers()!Norm(S_c);
                printf "sums of roots of unity for c = %o: %o\n", c, Nm_s ne 0 select PrimeDivisors(Nm_s) else [0];
                n_c := Minimum([n : n in [1..Nm_s] | n * OK subset S_c * OK]);
                printf "n_c = %o\n", n_c;
                lcm_n_c := LCM(lcm_n_c, n_c);
                //printf "sums of roots of unity for c = %o: %o\n", c, Nm_s ne 0 select Factorization(ideal<OK|Nm_s>) else [0] where Nm_s := (OK!(sum));
                //printf "sums of roots of unity for c = %o: %o\n", c, MinimalPolynomial(Nm_s) where Nm_s := (OK!(sum));
            end for;
            printf "Lcm(n_c) = %o\n\n", lcm_n_c;
            end for;
        else
            //printf "J_0(%o)^new has no rank 0 quotient.\n", M;
        end if;
    end for;
    printf "\n";
end for;

/*
file := Open("xqyq.sage", "w");
for N in minimal_levels do
    J := J0Nstar(N);
    g := Dimension(J);
    //r, r0quotient := AnalyticRanks(J : find_rank_0_quotient := false);
    //printf "N = %o: g_0(N)^* = %o, r = %o (has rank 0 quot: %o; r < g: %o)\n", N, g, r, r0quotient, r lt g;
    label := AnalyticRanks(J : find_rank_0_quotient := false);
    if Type(label) eq RngIntElt then
        printf "N = %o has no rank 0 quotient.\n", N;
        continue;
    end if;
    NE := -1;
    for i -> NElabel in label do 
        if NElabel[1] gt NE then 
            NE := NElabel[1];
            E := EllipticCurve(NElabel[2]);
        end if;
    end for;
    if NE eq -1 then 
        printf "no E found!\n";
    else
        printf "N = %o: using E = %o\n", N, CremonaReference(E);
        fprintf file, "E = EllipticCurve(\"%o\")\n", CremonaReference(E);
        fprintf file, "xq, yq = E.modular_parametrization().power_series(prec = 5*(%o + 1))\n", NE;
        fprintf file, "print(\"[* {}, '{}', {}, {}, {} *], \".format(%o, '%o', xq, yq, %o))\n\n", NE, CremonaReference(E), N;
    end if;
end for;

// in Geany: replace
//  \+ O\(q\^\d+\)
// by the empty string
// and remove the last ,
// replace ' by "

// sanity check:
assert forall{N : N in minimal_levels | N notin g0Nstarnull};

exceptional_levels := [N : N in [2..100000] | IsExceptional(N)];
for wtN in exceptional_levels do 
    if not exists{N : N in minimal_levels | IsDivisibleBy(wtN, N) and IsSquare(ExactQuotient(wtN, N))} then 
        if wtN in g0Nstarnull then
            continue;
        end if;
        printf "N = %o = %o is exceptional, but not a square times a minimal exceptional level; g_0(N)^* = %o\n", wtN, Factorization(wtN), GenusX0NStar(wtN);
    end if;
end for;*/