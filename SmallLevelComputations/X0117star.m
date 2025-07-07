load "../QuadraticPoints/models_and_maps.m";
load "../J0wplusminus.m";
load "../gonal_maps.m";
load "../Coleman/coleman.m";
load "../modelsX0Nstar.m";
SetDebugOnError(true);


N := 117;

function classify_points_on_X0N_wd(N, wds)
    time X, ws, pairs, NB, cusp := eqs_quos(N, [wds]);
    printf "X_0(%o)/%o: %o\n", N, wds, pairs[1][1];

    time j := jmap(X,N);
    printf "computed j: X_0(%o) -> X(1).\n", N;

    X0Nw, quotient := Explode(pairs[1]);
    pts := RationalPoints(X0Nw : Bound := 1000);
    printf "%o small points: %o\n", #pts, pts;

    function printj(pt)
        x := pt[1];
        y := pt[2];
        if y eq 0 then
            return "cusp";
        else
            assert y eq 1;
            try
                x := Rationals()!x;
                return Sprint(x);
            catch e 
                return Sprint((MinimalPolynomial(x)));
            end try;
        end if;
    end function;

    eqnsquo := DefiningEquations(quotient);
    eqnsj := DefiningEquations(j);

    for pt in pts do
        pts := PointsOverSplittingField(pt@@quotient); // only [1]: all points correspond to isogenous curves
        printf "pt = %o: \n j's = %o\n", pt, {printj(j(ptj)) : ptj in pts};
    end for;

    return X0Nw;
end function;


function reduce(P, Xp)
    return Divisor(Xp!ChangeUniverse(Eltseq(P), BaseField(Xp)));
end function;

procedure finite_index_subgroup(X, pts, r)
    bp := pts[1];
    pts := pts[2..#pts];
    if #pts lt r then 
        printf "#pts-1 < r\n";
        return;
    end if;
    for p in PrimesUpTo(50) do 
        try
            Xp := ChangeRing(X, GF(p));
            PicXp, phi, psi := ClassGroup(Xp);
            JFp := TorsionSubgroup(PicXp);
            divsp := [psi(reduce(P,Xp) - reduce(bp,Xp)) : P in pts];
            mat := Matrix(Integers(), #divsp, #Generators(JFp), [Eltseq(JFp ! D) : D in divsp]);
            for ell in PrimesUpTo(40) do
                if Rank(ChangeRing(mat, GF(ell))) eq r then
                    printf "p = %o: ell = %o\n", p, ell;
                    //break p;
                    return;
                end if;
            end for;
        catch e;
            //print e;
        end try;
    end for;
end procedure;

procedure rational_points(N, d : X := 0)
    if X cmpeq 0 then
        wds := [[1, d]];
        time X, ws, pairs, NB, cusp := eqs_quos(N, wds);
        printf "X_0(%o)/w_%o: %o\n", N, d, pairs[1][1];
        X, quotient := Explode(pairs[1]);
    end if;

    pts := PointSearch(X, 100);
    printf "%o small points: %o\n", #pts, pts;

    // rank of X0Nw
    J0N := JZero(N);
    J0Nwd := ConnectedKernel(1 - AtkinLehnerOperator(J0N, d)); // sign +1
    g := Dimension(J0Nwd);
    printf "g = %o\n", g;
    r := AnalyticRanks(J0Nwd);
    printf "r = %o\n", r;

    assert r lt g;
    printf "TODO: differences of %o points generate rank %o subgroup of J_0(N)^w(Q).\n", #pts-1, r;
    finite_index_subgroup(X, pts, r);

    // find a plane model
    try
        Qs := [];
        denoms := [];
        mindegree := 100;
        for i -> bp in pts do
            //model1, m := ProjectionFromNonsingularPoint(X, bp);
            //plane_model := ProjectionFromNonsingularPoint(model1, m(pts[i lt #pts select i + 1 else 1]));
            bp2 := bp;
            plane_model := X;
            while Dimension(AmbientSpace(plane_model)) gt 2 do
                plane_model := ProjectionFromNonsingularPoint(plane_model, bp2);
                bp2 := PointSearch(plane_model, 100)[1];
            end while;
            image_curve_non_monic_eq_xy := Evaluate(DefiningEquation(plane_model), [x, y, 1]);
            image_curve_non_monic_eq_xy := image_curve_non_monic_eq_xy / LeadingCoefficient(LeadingCoefficient(image_curve_non_monic_eq_xy));
            image_curve_non_monic_eq_xy_lc := LeadingCoefficient(image_curve_non_monic_eq_xy);
            Q := Numerator(Evaluate(image_curve_non_monic_eq_xy, y / image_curve_non_monic_eq_xy_lc));
            //Q := Evaluate(Q, [10*x,y]);
            //print Q, "\n";
            Append(~Qs, Q);
            Append(~denoms, Lcm([Lcm([Denominator(Coefficient(Qx, j)) : j in [0..Degree(Qx)]]) : Qx in Coefficients(Q)]));
            mindegree := Min(mindegree, Degree(Q));
        end for;

        print Qs, denoms;
    catch e 
        print e;
    end try;
    min_denom, min_i := Minimum(denoms);
    i := min_i;
    Q := Qs[i];
    tempQ := Qs[1];
    Qs[1] := Q;
    Qs[i] := tempQ;

    for i -> Q in Qs do
        printf "\ntrying model #%o: %o\n", i, Q;
        // Chabauty
        for p in PrimesInInterval(2, 50) do
            try 
                time data := coleman_data(Q, p, 25);
                time L,v := effective_chabauty(data : bound := 1000, e := 50); 
                QXYZ<X,Y,Z> := PolynomialRing(Rationals(), 3);
                mQ := &+[&+[Coefficient(Qx, j) * X^j : j in [0..Degree(Qx)]] * Y^(i-1) : i -> Qx in Coefficients(Q)];
                mQ := Homogenization(mQ, Z);
                rpts := RationalPoints(Curve(ProjectiveSpace(QXYZ), mQ) : Bound := 1000);
                printf "L = %o (%o pts)\nratpts = %o (%o pts)\n", L, #L, rpts, #rpts;
                if #L eq #rpts then
                    printf "found all Q-points!\n";
                    return;
                else
                    printf "have to exclude %o residue discs.\n", #L - #rpts;
                end if;
            catch e 
                printf "p = %o fails: %o.\n", p, e;
            end try;    
        end for;
    //    break;
    end for;
end procedure;

SetLogFile("X0117.log");
printf "X_0(117)^*:\n";
wds := [d : d in Divisors(N) | Gcd(d, ExactQuotient(N,d)) eq 1];
X := classify_points_on_X0N_wd(N, wds);
