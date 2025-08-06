load "../QuadraticPoints/models_and_maps.m";
load "../J0wplusminus.m";
load "../Coleman/auxpolys.m";
load "../Coleman/coho.m";
load "../Coleman/froblift.m";
load "../Coleman/reductions.m";
load "../Coleman/singleintegrals.m";
load "../Coleman/applications.m";
load "../modelsX0Nstar.m";
SetDebugOnError(true);

N := 396;

SetLogFile("X0396star" * index * ".log");
index := StringToInteger(index);

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
                    // compute linearly independent divisors[1..r]
                    // splitting_indices = [[0..1..0] : i in [1..r]]
                    return;
                end if;
            end for;
        catch e;
            //print e;
        end try;
    end for;
end procedure;

procedure rational_points(N)
    X := XZeroNstar(N);

    pts := PointSearch(X, 100);
    printf "%o small points: %o\n", #pts, pts;

    // rank of X0Nw
    J0Nstar := J0wplusminus([p^Valuation(N,p) : p in PrimeDivisors(N)], [Integers()|]);
    g := Dimension(J0Nstar);
    printf "g = %o\n", g;
    r := AnalyticRanks(J0Nstar);
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
                plane_model, quot := ProjectionFromNonsingularPoint(plane_model, bp2);
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

        //print Qs, denoms;
    catch e 
        print e;
    end try;
    min_denom, min_i := Minimum(denoms);
    i := min_i;
    Q := Qs[i];
    tempQ := Qs[1];
    Qs[1] := Q;
    Qs[i] := tempQ;

    for i -> Q in Qs[index..index] do
        printf "\ntrying model #%o: %o\n", i, Q;
        // Chabauty
        for p in PrimesInInterval(2, 50) do
            try 
                time data := coleman_data(Q, p, 75);
                time L,v := effective_chabauty(data : bound := 1000, e := 50);
                printf "found %o residue disks.\n", #L;
                /*QXYZ<X,Y,Z> := PolynomialRing(Rationals(), 3);
                mQ := &+[&+[Coefficient(Qx, j) * X^j : j in [0..Degree(Qx)]] * Y^(i-1) : i -> Qx in Coefficients(Q)];
                mQ := Homogenization(mQ, Z);
                rpts := RationalPoints(Curve(ProjectiveSpace(QXYZ), mQ) : Bound := 1000);*/
                rpts := Q_points(data, 1000);
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


rational_points(N);
