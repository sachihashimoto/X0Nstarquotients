//https://warwick.ac.uk/fac/cross_fac/complexity/people/students/dtc/students2013/klaise/janis_klaise_ug_report.pdf

omega := function(N)
    return #[p : p in PrimeDivisors(N)];
end function;


function maxDisc(powerof2)
        if powerof2 eq 0 then 
            return 163; //this is the largest |D_O| where O class number 1
        elif powerof2 eq 1 then 
            return 427; //largest |D_O|, O class number 2
        elif powerof2 eq 2 then 
            return 1555; //largest |D_O|,  O class number 4
        elif powerof2 eq 3 then 
            return 7987;
        elif powerof2 eq 4 then 
            return 35275;
        else
            error "not implemented";
        end if;
end function;


//(1) For any N in our short list, find all pairs (D_K,c) such that 
// the order of discriminant D_O = D_K c^2 has Picard group elementary 2-torsion of order dividing 2^{\omega(N)} 

function findOrders(N)
    omegaN := omega(N);
    orders := [];
    for Dplus in [3.. maxDisc(omegaN)] do 
    //iterate through all possible discriminants of orders
        D := -Dplus;
        K := QuadraticField(D);
        OK := Integers(K);
        DK := Discriminant(K);
        
        if not IsDivisibleBy(D, Discriminant(K)) then
            continue;
        end if;
        O := sub<OK | Round(Sqrt(ExactQuotient(D, DK)))>;
        c := Conductor(O);
        assert DK*c^2 eq D;
        assert Discriminant(O) eq D;

        ClO, mClO := PicardGroup(O);
        if IsDivisibleBy(2^omegaN, #ClO) then
            // check that it's an *elementary* 2-group
            if &and[e+e eq Identity(ClO) : e in ClO ] then
                Append(~orders, [DK, c]);
            end if;
        end if;
    end for;
    return orders;
end function;


function FormClass(I)
//given an ideal I in an order O of discriminant D
//return a form in the form class group representing I
    g1, g2 := Explode(Basis(I));
    O := Order(I);
    D := Discriminant(O);
    f := Conductor(O);
    K := QuadraticField(D);
    _<z> := PolynomialRing(K);
    DK := Discriminant(K);
    rts := Roots(z^2 - DK);
    sqrtDK := rts[1][1];
    
    //want Norm(xg1 - yg2)/Norm(I)
    NI := Norm(I);
    _<x> := PolynomialRing(K);
    g := x*g1- g2;
    gconj :=x*Conjugate(g1) -Conjugate(g2);
    Nmg :=  g*gconj;
    e1,e2,e3 := Explode(Coefficients(Nmg/NI));
    den := LCM([Denominator(e1), Denominator(e2), Denominator(e3)]);
    Q:=QuadraticForms(D);
    cl := Q![Integers()!(e3*den),Integers()!(e2*den),Integers()!(e1*den)];
    return ReducedForm(cl);
end function;

//given N, D_K and c fixed, list (in increasing order of the prime factors) for each p^k || N 
//the admissible prime ideals of norm p^k 
//(grouped by conjugate pairs, with an arbitrary numbering for each) << not currently grouped by conjugate pairs

function AdmissiblePrimeIdealsConditionA(N, pair)
    DK, c := Explode(pair);
    primesN := PrimeFactors(N);
    coprime_to_c := [p : p in primesN | not(IsDivisibleBy(c, p))];
    D := DK *c^2;
    K := QuadraticField(DK);
    OK := Integers(K);
    O := sub<OK | Round(Sqrt(ExactQuotient(D, DK)))>;
    I := AssociativeArray();
    Omax := MaximalOrder(K); //magma can only factor over the maximal order
    //factor over maximal order, then restrict back to O
    if coprime_to_c eq [] then
        return true, [];
    end if;
    for p in coprime_to_c do
        good := IsSplit(p, OK) or (IsRamified(p, OK) and Valuation(N, p) eq 1); //check condition (a)
        if not good then
            return false, [];
        end if;
         v := Valuation(N, p); //need to see how many copies we want
         //for each copy we can choose either ideal the plus or minus
         fact := [];
         for fac in Factorization(p*Omax) do
             if fac[2] eq 2 then
               // fact := fact cat [fac[1],fac[1]];
               fact := fact cat [fac[1]];
             else
                 assert fac[2] eq 1;
                 fact := fact cat [fac[1]];
             end if;
         end for;
         //assert #fact eq 2;
         Ip := [fac meet O : fac in fact]; //now to restrict to O-ideal
         //assert #Ip eq 2;
         I[p] := [{* Ip[i]^^v *} : i in [1 .. #Ip]]; //1 or 2 options depending on ramified or split
    end for;
    flattenI := [I[p] : p in Keys(I)];

    //flattenI is a sequence of elements of the cartesian product
    //need to make sequence of multisets
 
    return good, flattenI;
    //return a list of lists
    //each list contains all prime ideals of norm p^k
    //a prime ideal of norm p^k is given as a multiset of prime ideals of norm p
end function;




function AdmissiblePrimeIdealsConditionB(N, pair)
    DK, c := Explode(pair);
    primesN := PrimeFactors(N);
    coprime_to_c := [p : p in primesN | not(IsDivisibleBy(c, p))];
    D := DK *c^2;
    K := QuadraticField(DK);
    OK := Integers(K);
    O := sub<OK | Round(Sqrt(ExactQuotient(D, DK)))>;
    I := AssociativeArray();
    Omax := MaximalOrder(K); 
    conditionb := true;
    //magma can only factor over the maximal order
    //factor over maximal order, then restrict back to O
    m := GCD(N, c);
    primes_m := PrimeDivisors(m);
    for p in primes_m do
        k := Valuation(N, p);
        lambdaset := [];
        if DK mod 4 eq 0 then
            Dprime := c^2 * DK/4;
            alpha := Sqrt(O!(Dprime));
            for lambda in [0..p^k-1] do
                if (lambda^2 - Integers()!Dprime) mod p^k eq 0 and Valuation(Dprime-lambda^2, p) eq k then
                    Append(~lambdaset,lambda);
                end if;
            end for;
        else
            assert DK mod 4 eq 1;
            alpha := (O!c+Sqrt(O!(c^2*(DK))))/2;
            for lambda in [0..p^k-1] do
                //solve equation, valuation
                if Valuation(lambda^2 + lambda*c - c^2 * ExactQuotient(DK-1,4), p) eq k then
                    Append(~lambdaset,lambda);
                end if;
            end for;

        end if;
        if #lambdaset eq 0 then
            conditionb := false;
            return conditionb, [];
        end if; 
        if p in Keys(I) then
            for i in [1..#lambdaset] do
                Append(~I[p], ideal<O |O!p^k, alpha + (O!lambdaset[i])> );
            end for;
        else
        I[p] := [ideal<O |O!p^k, alpha + (O!lambdaset[i]) >: i in [1..#lambdaset]]; //I[p] contains the admissible ideals of norm p^k
        end if;
    end for;
    if #primes_m eq 0 then
        conditionb := true;
        return conditionb, []; //condition b is true, but there are no admissible ideals
    end if;
    flattenI := [ I[p]: p in Keys(I)]; 
    //return a list of lists, each list corresponds to a prime p, and contains admissible ideals of norm p^k
    fixtypeI := [ [{* ip *} : ip in i ] : i in flattenI];
    
    return conditionb, fixtypeI;

end function;


function AdmissiblePrimeIdeals(N, pair : verbose := true)
            b, A := AdmissiblePrimeIdealsConditionA(N,pair);
            if b then
                b, B := AdmissiblePrimeIdealsConditionB(N, pair);
                    if b then
                        if verbose then
                            printf "Pair %o does satisfy both conditions.\n", pair;
                        end if;
                        //now pair passes condition A and condition B
                        //have the sequences A and B of admissible ideals
                        ideals := []; // make a list of lists, indexed by primes dividing N
                        //each list contains all admissible ideals of norm p^k, where k = valp(N)
                        primesN := PrimeFactors(N);
                        DK, c := Explode(pair);
                        m := GCD(N, c);
                        primesm := PrimeFactors(m);
                        for p in primesN do
                            k := Valuation(N, p);
                            if p in primesm then
                                //this is in list B
                                listb := [b: b in B |Norm(MultisetToSequence(b[1])[1]) eq p^k ]; //should be unique
                                assert #listb eq 1;
                                ideals_norm_pk := listb[1];
                            else
                                //this is in list A
                                lista := [a: a in A | Norm(MultisetToSequence(a[1])[1]) eq p];      //there are multiple options
                                assert #lista eq 1;
                                ideals_norm_pk := lista[1];
                            end if;
                            if verbose then
                                printf "The admissible prime ideals of norm %o^%o are:\n",p, k;
                                print ideals_norm_pk;
                            end if;
                            Append(~ideals, ideals_norm_pk);
                        end for;
                        return true, ideals;
                    else
                        if verbose then
                            printf "Pair %o does not satisfy condition (b) for %o.\n", pair, N;
                        end if;
                        return false, [];
                    end if;
        else
            if verbose then
                printf "Pair %o does not satisfy condition (a) (the Heegner hypothesis) for %o.\n", pair, N;
            end if;
            
            return false, [];    
        end if;

end function;


//function for dealing with multisets instead of ideals
function MultisetToIdeal(etap)
    //take multiset extract elements, multiply together
    etap := MultisetToSequence(etap);
    I := &*[elt : elt in etap];
    return I;
end function;

function FindComplexConjugate(etap)
//let etap be a multiset denoting an admissible ideal of norm p^k
//return the complex conjugate
    etap := MultisetToSequence(etap);
    conj := [Conjugate(elt) : elt in etap];
    return Multiset(conj);
end function;

function NormOfMultiset(etap)
     etap := MultisetToSequence(etap);
     I := &*etap;
     return Norm(I);
end function;

function FormClassMultiset(etap)
    etap := MultisetToSequence(etap);
    I := &*etap;
    return FormClass(I);
end function;


//(2) //for each of those, and representatives of all elements of the class group, give its representatives
//(3) //Now, a Heegner point is represented as follows : 
//((D_K,c),(eta_p)_p,\mathfrak{a}) where \mathfrak{a} represents a class in the class group.

function HeegnerPointsOrder(pair, N : verbose := true)
    HPs := [**];
        b, A := AdmissiblePrimeIdeals(N,pair : verbose := verbose);
        if b then
            //indexed by p dividing N
            etas := CartesianProduct(A);
            forms := ReducedForms(pair[1]*pair[2]^2);
            for eta in etas do
                //each eta correspons to a choice of ideal of norm p^k for each p
                //assign the form class
                for f in forms do 
                    hp := <pair, eta, f>;
                    Append(~HPs, hp);
                end for;
            end for;
        end if;
    return HPs;
end function;


function HeegnerPoints(N: verbose:= false)
    pairs := findOrders(N);
    HPs := [**];
    for pair in pairs do 
        b, A := AdmissiblePrimeIdeals(N,pair : verbose := verbose);
        if b then
            //indexed by p dividing N
            etas := CartesianProduct(A);
            forms := ReducedForms(pair[1]*pair[2]^2);
            for eta in etas do
                //each eta correspons to a choice of ideal of norm p^k for each p
                //assign a form class
                for f in forms do 
                    hp := <pair, eta, f>;
                    Append(~HPs, hp);
                end for;
            end for;
        end if;
    end for;
    return HPs;
end function;


//w_Q simply switches the admissible prime ideals between conjugate pairs for p|Q, 
///does not nothing for the others, and multiplies by product of classes of the [eta_p] (p|Q) in the class group.


function HallDivisors(N) 
    return [q : q in Divisors(N) | Gcd(q, ExactQuotient(N, q)) eq 1];
end function;

function actAtkinLehner(hp, Q, N)
//for a single heegner point, return w_Q^(N) of hp 
    pair := hp[1];
    eta := hp[2];
    a := hp[3];
    etaprime := <>; //this is the new eta
    etaQ := [];
    for etap in eta do
        pk := NormOfMultiset(etap);
        if Integers()!Q mod Integers()!pk eq 0 then 
            etapconj := FindComplexConjugate(etap);
            Append(~etaprime, etapconj);
            Append(~etaQ, etap);
        else
            Append(~etaprime, etap);
        end if;
    end for;

    if #PrimeFactors(Q) eq 0 then 
        return <pair, eta, a>;
    else
        //need to change the formclass
        fQ := &*[FormClassMultiset(etap) : etap in etaQ];
    end if;
    return <pair, etaprime, ReducedForm(a*fQ)>;
    
end function;



function simpleForm(hp)
    //give the heegner point but with unfactored form in the middle
    eta := hp[2];
    I := &*[MultisetToIdeal(etap) : etap in eta];
    return <hp[1], I, hp[3]>;
end function;

function checkdone(hp,done)
    if #done eq 0 then
        return false;
    end if;
    for elt in done do
        areequal := true;
        for i in [1..3] do
            if (hp[i] cmpne elt[i]) then
                areequal := false;
            end if;
        end for;
        if areequal then
            return true;
        end if;
    end for;
    return false;
end function;


//the complex conjugate switches between conjugate pairs for all p|N, and multiplies by the classes of the [eta_p] (all p).
//An element of Gal(H/K) corresponding to a class just multiplies inside the class group, does not touch the (eta_p)_p part.

function tau(hp) 
    //this function should return hp acted on by complex conjugation
    etaconj := <FindComplexConjugate(etap) : etap in hp[2]>;
    return <hp[1], etaconj, hp[3]^(-1) >;
end function;

function galoisOrbit(hp)
    //this is the action by all elements of Pic(O) as well as complex conj.
    pair := hp[1];
    DK, c := Explode(pair);
    D := DK *c^2;
    ClD := ReducedForms(D);  //list of forms in form class group 
    done := [];
    galorbit := [];
    for f in ClD do 
        hpf :=  <hp[1],hp[2], ReducedForm(hp[3] *  f) >;
        b := checkdone(simpleForm(hpf), done);
        if b then
            continue;
        end if;
        Append(~done, simpleForm(hpf));
        Append(~galorbit, hpf);
    end for;
    
    hpconj := tau(hp);
    
    for f in ClD do 
        hpf := <hpconj[1],hpconj[2], ReducedForm(hpconj[3]*f)>;
        b := checkdone(simpleForm(hpf), done);
        if b then
            continue;
        end if;
        Append(~done, simpleForm(hpf));
        Append(~galorbit, hpf);
    end for;
    
    return galorbit;
end function;



function allGaloisOrbits(hps)
    //return all galois orbits
    done := [* *];
    equivalence_classes := [**];
    for i->hp in hps do
        b := checkdone(simpleForm(hp), done);
        if b then 
            continue;
        end if;
        orbit := galoisOrbit(hp);
        for elt in orbit do
            Append(~done, simpleForm(elt));
        end for;
        Append(~equivalence_classes, orbit);
            
    end for;
    
    return equivalence_classes;

end function;


function equivalenceClasses(hps, N)
    //return equivalence classes of heegner points under Atkin--Lehner
    done := [* *];
    equivalence_classes := [**];
    for i->hp in hps do
        equiv_hp := [**] ;
        b := checkdone(simpleForm(hp), done);
        if b then 
            continue;
        end if;
        Append(~equiv_hp, hp);
        Append(~done, simpleForm(hp));
        for Q in HallDivisors(N) do
            AL := actAtkinLehner(hp, Q, N);
            b := checkdone(simpleForm(AL), done);
            if b then 
                continue;
            end if;
            Append(~equiv_hp,AL);
            Append(~done, simpleForm(AL));
        end for;
        Append(~equivalence_classes, equiv_hp);
    end for;
    
    return equivalence_classes;
end function;
//need to check that the Galois orbits are contained in a single AL orbit

function isContainedIn(L1, L2);
    //check containment of L1 in L2.
    for l1 in L1 do
        are_equal := false;
        for l2 in L2 do
            if l1 cmpeq l2 then 
                are_equal := true;
                continue;
            end if;
        end for;
        if not are_equal then
            return false;
        end if;
    end for;
    return true;
end function;



function galoisALcompatibleHps(hps, N) 
// this only checks if they are all compatible
//we want to check which ones are compatible, and throw out all hps associated with non-compatible orders
// if i'm correct, this returns the rational points on X0(N)*?
    E_al := equivalenceClasses(hps, N);
    E_g := allGaloisOrbits(hps);
    badorders := [];
    for orbit in E_g do
        pair :=orbit[1][1];
        isOK := &or[isContainedIn(orbit, eqclass)  : eqclass in E_al]; 
        if not isOK then
            ind := Index(E_g, orbit);
            Remove(~E_g, ind);
        end if;
    end for;
    pts := [];
    for orbit in E_g do
        for eqclass in E_al do
            if isContainedIn(orbit, eqclass) then 
                Append(~pts, eqclass);
                ind := Index(E_al, eqclass);
                Remove(~E_al, ind);
            end if;
        end for;
    end for;
    return pts;
end function;

procedure myprint(pts)
    for pt in pts do
        a, b := Explode(pt[1][1]);
        print(a*b^2); //CM order disc
    end for;
end procedure;


function all_discs(pts)
    discs := [];
    for pt in pts do
        a, b := Explode(pt[1][1]);
        Append(~discs, (a*b^2)); //CM order disc
    end for;
    return discs;
end function;


