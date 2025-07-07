load "../heegner.m";

for N in [Integers()| 52/4, 68/4, 164/4, 148/4, 260/4, 212/4] do
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);
    Ds := all_discs(allpts);
    assert -16 in Ds;
end for;


for N in [Integers()| 63/9, 171/9, 279/9] do
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);
    Ds := all_discs(allpts);
    assert -27 in Ds;
end for;

for N in [Integers()|84/4, 156/4, 228/4, 444/4] do
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);
    Ds := all_discs(allpts);
    assert -48 in Ds;
end for;

for N in [Integers()|98/2, 294/2] do
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);
    Ds := all_discs(allpts);
    assert -12 in Ds;
end for;

for N in [Integers()|294/2] do
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);
    Ds := all_discs(allpts);
    assert -48 in Ds;
end for;


for N in [Integers()|242/2] do
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);
    Ds := all_discs(allpts);
    assert -28 in Ds;
end for;

for N in [Integers()|308/4] do
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);
    Ds := all_discs(allpts);
    assert -112 in Ds;
end for;





