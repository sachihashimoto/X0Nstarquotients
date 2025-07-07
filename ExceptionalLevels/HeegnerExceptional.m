SetLogFile("HeegnerPointsExceptionalLevels.log");
load "../heegner.m";

for N in [40,48,72,80,88,96,100,108,112,120,135,144, 147, 162, 176,180,184, 196,200,216,224,225,240,396] do
	print "Level";
	print N;
	hps := HeegnerPoints(N);
	allpts := galoisALcompatibleHps(hps, N);

	for pt in allpts do
		a, b := Explode(pt[1][1]);
    	print(a*b^2); //CM order disc
    end for;
end for;
