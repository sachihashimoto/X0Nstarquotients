function J0wplusminus(wplus, wminus)
    J0Nwplus := JZero(&*wplus * &*wminus);
    for pk in wplus do
        J0Nwplus := ConnectedKernel(1 - AtkinLehnerOperator(J0Nwplus, pk)); // sign +1
    end for;
    J0Nwplusminus := J0Nwplus;
    for qk in wminus do
        J0Nwplusminus := ConnectedKernel(1 + AtkinLehnerOperator(J0Nwplusminus, qk)); // sign -1
    end for;
    return J0Nwplusminus;
end function;

function GenusX0NStar(N)
    J := JZero(N);
    seq_al := [qs[1]^qs[2] : qs in Factorization(N)];

    M:=CuspidalSubspace(ModularForms(N,2)); d:=2; dim:=Dimension(M);
    if dim eq 0 then return 0; end if;
    psi:=N*&*[1+1/f[1] : f in Factorization(N)]; prec:=Max(200,Ceiling(N*4*2/12));
    B:=qExpansionBasis(M,prec);
    S<q>:=Parent(B[1]); HYPER:=false;
    if #seq_al ne 0 then
     AL:=[AtkinLehnerOperator(M,al) : al in seq_al]; J,T:=Diagonalization(AL);
     B:=ChangeRing(Matrix([B]),PowerSeriesRing(Rationals()));
     B:=Eltseq(ChangeRing(T,Parent(B[1][1]))*Transpose(B));
     B:=[B[i] : i in [1..#B] | &and[a[i][i] eq +1 : a in J]];
    end if;
        return #B;

 //   X := X0NQuotient(N, seq_al);
 //   return Type(X) eq Prj select 0 else Genus(X);
    //return IsIsomorphic(X, CProjectiveSpace(Rationals(),1)) select 0 else Genus(X);
    //return Dimension(Image(&*[1 + AtkinLehnerOperator(J, i) : i in seq_al]));
    //return Dimension(J0wplusminus([qs[1]^qs[2] : qs in Factorization(N)], [Integers()|]));
end function;

function J0Nstar(N)
    J := JZero(N);
    return Image(&*[1 + AtkinLehnerOperator(J, i) : i in [qs[1]^qs[2] : qs in Factorization(N)]]);
end function;

function HasNewRank0Quotient(A)
    if Dimension(A) eq 0 then
        return false, [], [];
    end if;
    N := Level(A);
    As := Decomposition(A);
    rank0s := [];
    JMs := [];
    rank0 := false;
    for B in As do
        if Level(ModularSymbols(B)[1]) ne Level(ModularSymbols(A)[1]) then 
            continue;
        end if;
        L_is_zero_at_1 := IsZeroAt(LSeries(B), 1); //Abs(Evaluate(LSeries(B), 1)) lt 10^-7;
        if not L_is_zero_at_1 then
            rank0 := true;
        end if;
        Append(~rank0s, <not L_is_zero_at_1, Level(ModularSymbols(B)[1]), Dimension(B)>);
        Append(~JMs, B);
    end for;
    return rank0, rank0s, JMs;
end function;

function HasRank0Quotient(A)
    if Dimension(A) eq 0 then
        return false, [];
    end if;
    N := Level(A);
    As := Decomposition(A);
    rank0s := [];
    rank0 := false;
    for B in As do
        L_is_zero_at_1 := IsZeroAt(LSeries(B), 1); //Abs(Evaluate(LSeries(B), 1)) lt 10^-7;
        if not L_is_zero_at_1 then
            rank0 := true;
        end if;
        Append(~rank0s, <not L_is_zero_at_1, Level(ModularSymbols(B)[1]), Dimension(B)>);
    end for;
    return rank0, rank0s;
end function;


function AnalyticRanks(A : find_rank_0_quotient := true)
    if Dimension(A) eq 0 then
        //printf "0";
        return 0, false;
    end if;
    N := Level(A);
    As := Decomposition(A);
    if not find_rank_0_quotient then
        //printf "Ranks: ";
    end if;
    found_rank_0_quotient := false;
    r := 0;
    Es := [];
    for B in As do
        //printf "%o, ", IsZeroAt(LSeries(A), 1) select Dimension(A) else 0;
        // Abs(Evaluate(LSeries(A), 1)) lt 10^-7 is faster, but not provably correct
        L_is_zero_at_1 := Abs(Evaluate(LSeries(B), 1)) lt 10^-7;
        //L_is_zero_at_1 := IsZeroAt(LSeries(A), 1);
        if not L_is_zero_at_1 then
            found_rank_0_quotient := true;
            if Dimension(B) eq 1 then
                E := EllipticCurve(B);
                NE := Conductor(E);
                //good := not IsSquare(ExactQuotient(N, NE));
                //printf "%o (%o), ", CremonaReference(E), good;
         //       printf "%o, ", CremonaReference(E);
                //if good then
                    Append(~Es, <NE, CremonaReference(E)>);
                //end if;
                //return CremonaReference(EllipticCurve(A));
                //print CremonaReference(E), MordellWeilGroup(E);
            end if;
        end if;
        //print L_is_zero_at_1 select 1 else 0;
        r +:= L_is_zero_at_1 select Dimension(B) else 0;
        if not find_rank_0_quotient then
            printf "%o (N = %o), ", L_is_zero_at_1 select Dimension(B) else 0, Level(ModularSymbols(B)[1]);
        end if;
    end for;
    /*if not find_rank_0_quotient then
        return r, found_rank_0_quotient;
        //printf "(r = %o: rank 0 quotient: %o)", r, found_rank_0_quotient;
    end if;*/
    //printf "r = %o\n", r;
    return r;//return Es;
end function;