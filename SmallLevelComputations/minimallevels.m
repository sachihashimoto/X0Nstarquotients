function FilterExceptionalLevelsA(nprimes, p, exp)
    //return all nprime that are not exceptional for p^exp
    exceptional1 := [72, 100, 108, 144, 180, 196, 200, 216, 225, 324, 441, 450, 500, 1125, 1225, 1372];
    exceptional2 := [ 40, 48, 80, 88, 96, 99, 112, 120, 135, 147, 162, 169, 176, 184, 224, 240, 250, 297, 368, 405, 486, 1029, 1250];
    E := exceptional1 cat exceptional2;
    return [m : m in nprimes | p^exp * m notin E];
end function;


function FilterExceptionalLevelsB(levels)
    //return all n that are not exceptional
    exceptional1 := [72, 100, 108, 144, 180, 196, 200, 216, 225, 324, 441, 450, 500, 1125, 1225, 1372];
    exceptional2 := [ 40, 48, 80, 88, 96, 99, 112, 120, 135, 147, 162, 169, 176, 184, 224, 240, 250, 297, 368, 405, 486, 1029, 1250];
    E := exceptional1 cat exceptional2;
    return [m : m in levels | m[1] notin E];
end function;

function CreateNprime(B, p,exp,sqfreelist)
    //create all nonexceptional levels N' with factors in of PrimesUpTo(B) with p the powerful prime
    bd := Floor(B/p^exp);
    nprimes := [n : n in sqfreelist |n lt bd];
    nprimes := FilterExceptionalLevelsA(nprimes, p, exp); //return nprimes that are not exceptional for p
    return nprimes;
end function;

function SquareFreeList(B)
    //precompute a table of N' squarefree less than B/2^2
    sqfree := [];
    Bprime := Floor(B/2^2);
    for i in [1..Bprime] do
        if IsSquarefree(i) then
            Append(~sqfree, i);
        end if;
    end for;
    return sqfree;
end function;

function LevelsAboveN(N, B)
    //return all levels that differ by a square factor from N up to bound B
    primelist := PrimesUpTo(Floor(SquareRoot(B/N)));
    levelsabove := [p^2*N : p in primelist];
    return levelsabove;
end function;

function CheckGenus(N)
    //return genus of X0^*(N)
    m := Index(Gamma0(N));
    number_of_terms:=Floor(m*2/12);
    C := CuspForms(N);
    g := Dimension(C);
    if g eq 0 then return 0; end if;
    // start by simply diagonalising all AL involutions
    wds := [AtkinLehnerOperator(C, d) : d in Divisors(N) | d gt 1 and GCD(d, ExactQuotient(N,d)) eq 1];
    diag, D := Diagonalization(wds);
    i_s := [i : i in [1..Nrows(diag[1])] | forall{wd : wd in diag | wd[i,i] eq +1}];
    qbas := Basis(C, 2*number_of_terms);
    bas := [&+[D[i,j] * qbas[j] : j in [1..#qbas]] : i in i_s];
    return #bas;
end function;

function CachedGenus(N, A)
    if N notin Keys(A) then
            A[N] := CheckGenus(N);
    end if;
        return A[N], A;
end function;

B:= 1000;
//by Gonzalez, Bars, we know that g > 5 when N > 910
//This code produces all levels N which are minimal and less than B, along with their genus
A := AssociativeArray();
all_levels := [];
sqfreelist := SquareFreeList(B);
for p in PrimesUpTo(B) do
    for exp in [2,3] do
        Nprimelist := CreateNprime(B,p,exp,sqfreelist);
        levels := [[p^exp *Nprime, CachedGenus(p^exp*Nprime,A)] : Nprime in Nprimelist];
        for pair in levels do
            if pair[2] eq 0 then
                levelsabove := LevelsAboveN(pair[1],B);
                levelsabove := [[N,CachedGenus(N,A)]: N in levelsabove ];
                remaininglevels := [pair : pair in levelsabove | pair[2] eq 0 ];
                goodlevels := SetToSequence(Set(levelsabove) diff Set(remaininglevels));
                goodlevels := FilterExceptionalLevelsB(goodlevels);
                all_levels := all_levels cat goodlevels;
                while #remaininglevels gt 0 do
                    levelsabove := &cat [LevelsAboveN(R[1], B) : R in remaininglevels]; //go up
                    levelsabove := [[N,CachedGenus(N,A)]: N in levelsabove ]; //filter out
                    remaininglevels := [pair : pair in levelsabove | pair[2] eq 0 ];
                    goodlevels := SetToSequence(Set(levelsabove) diff Set(remaininglevels));                
                    goodlevels := FilterExceptionalLevelsB(goodlevels);
                    all_levels := all_levels cat goodlevels;
                end while;
            else
                Append(~all_levels, pair);
            end if;
        end for;
    end for;
end for;
all_levels := {x : x in all_levels | x[2] le 5};
print SetToSequence(all_levels);

