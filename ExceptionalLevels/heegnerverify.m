load "../heegner.m";

for N in [49,25] do
    print "Level";
    print N;
    hps := HeegnerPoints(N);
    allpts := galoisALcompatibleHps(hps, N);

    for pt in allpts do
        a, b := Explode(pt[1][1]);
        print(a*b^2); //CM order disc
    end for;
end for;
