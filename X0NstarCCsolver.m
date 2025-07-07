load "QuadraticPoints/models_and_maps.m";
load "J0wplusminus.m";
load "gonal_maps.m";
load "Coleman/coleman.m";
load "modelsX0Nstar.m";
SetDebugOnError(true);



function reduce(P, Xp)
    return Divisor(Xp!ChangeUniverse(Eltseq(P), BaseField(Xp)));
end function;

function finite_index_subgroup(X, pts, r)
    bp := pts[1];
    pts := pts[2..#pts];
    if #pts lt r then 
        error "#pts-1 < r\n";
    end if;
    for p in PrimesInInterval(5,15) do 
        try
            Xp := ChangeRing(X, GF(p));
            PicXp, phi, psi := ClassGroup(Xp);
            JFp := TorsionSubgroup(PicXp);
            divsp := [psi(reduce(P,Xp) - reduce(bp,Xp)) : P in pts];
            mat := Matrix(Integers(), #divsp, #Generators(JFp), [Eltseq(JFp ! D) : D in divsp]);
            for ell in PrimesUpTo(40) do
                if Rank(ChangeRing(mat, GF(ell))) eq r then
                    printf "finite index subgroup certified by: p = %o: ell = %o\n", p, ell;
                    return true;
                end if;
            end for;
        catch e;
            //print e;
        end try;
    end for;
    return false;
end function;

procedure rational_points(N , X )
    // rank of X0Nw
    J0Nstar := J0wplusminus([p^Valuation(N,p) : p in PrimeDivisors(N)], [Integers()|]);
    g := Dimension(J0Nstar);
    printf "g = %o\n", g;
    r := AnalyticRanks(J0Nstar);
    printf "r = %o\n", r;
    if r ge g then
        error "Chabauty condition not satisfied!";
    end if;

    if IsHyperelliptic(X) then
        pts := Points(SimplifiedModel(X) : Bound := 100);
        printf "checking: differences of %o points generate rank %o subgroup of J_0(N)^*(Q).\n", #pts-1, r;
        if not finite_index_subgroup(SimplifiedModel(X), pts, r) then 
            error "differences of small rat pts do not generate a finite index subgroup of J(Q).";
        end if;
    else
        pts := PointSearch(X, 100);
        printf "checking: differences of %o points generate rank %o subgroup of J_0(N)^*(Q).\n", #pts-1, r;
        if not finite_index_subgroup(X, pts, r) then 
            error "differences of small rat pts do not generate a finite index subgroup of J(Q).";
        end if;
    end if;
    printf "%o small points: %o\n", #pts, pts;




    // find a plane model
    try
        Qs := {};
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
            Include(~Qs, Q);
            Append(~denoms, Lcm([Lcm([Denominator(Coefficient(Qx, j)) : j in [0..Degree(Qx)]]) : Qx in Coefficients(Q)]));
            mindegree := Min(mindegree, Degree(Q));
        end for;

        //print Qs, denoms;
    catch e 
        print e;
    end try;
    printf "found %o plane models:\n%o\n", #Qs, Qs;

    for Q in Qs do
        printf "\ntrying model #\n%o", Q;
        // Chabauty
        for p in PrimesInInterval(4, 30) do
            try 
                time data := coleman_data(Q, p, 30);
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
N := 316;
SetLogFile("X0" * Sprint(N) * "star.log");
// N := StringToInteger(N);
print N;
//index := StringToInteger(index);
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
print X;

rational_points(N, X);
