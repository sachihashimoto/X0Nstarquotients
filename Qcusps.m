//short script to compute number of Q-cusps

function Qcusps(N)
	qcusps := 1;
	if N mod 4 eq 0 then
		qcusps +:= 1;
		if N mod 16 eq 0 and GCD(16, N div 16) eq 1 then
			qcusps +:=1;
		end if;
	end if;
	if N mod 9 eq 0 and GCD(9, N div 9) eq 1 then
		qcusps +:= 1;
		if N mod 4 eq 0 then
			qcusps +:=1;
			if N mod 144 eq 0 and GCD(144, N div 144) eq 1 then
				qcusps +:=1;
			end if;
		end if;
	end if;
	return qcusps;
end function;