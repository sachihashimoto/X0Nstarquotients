SetLogFile("HeegnerPointsRemaining.log");
load "../heegner.m";
for N in [ 63, 75, 81, 90, 98, 117, 121, 125, 126, 150, 153, 162, 171, 175, 189, 198, 207, 225, 234, 242, 245, 261, 270, 275, 279, 294, 306, 315, 342, 350, 378, 414, 495, 630 ] do
	print "Level";
	print N;
	hps := HeegnerPoints(N);
	allpts := galoisALcompatibleHps(hps, N);

	for pt in allpts do
		a, b := Explode(pt[1][1]);
    	print(a*b^2); //CM order disc
    end for;
end for;