SetLogFile("HeegnerPoints4.log");
load "../heegner.m";
for N in [ 220, 108, 68, 124, 48, 64, 84, 156, 140, 52, 80, 188, 132, 96, 76, 72, 100, 204, 276, 116, 180, 168, 112, 104, 380, 284, 420, 348, 312, 476, 152, 144, 240, 196, 136, 248, 252, 128, 164, 236, 300, 264, 172, 148, 280, 200, 176, 224, 260, 228, 308, 396, 316, 444, 376, 364, 216, 440, 208, 572, 212 ] do
	print "Level";
	print N;
	hps := HeegnerPoints(N);
	allpts := galoisALcompatibleHps(hps, N);

	for pt in allpts do
		a, b := Explode(pt[1][1]);
    	print(a*b^2); //CM order disc
    end for;
end for;

//Heegner points on small levels N where X0(N)* has genus g<6 and 4 | N 