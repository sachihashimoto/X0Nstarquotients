load "../J0wplusminus.m";
load "../gonal_maps.m";
load "../Coleman/coleman.m";
SetDebugOnError(true);

SetLogFile("X0Nstarhyperelliptic.log");

procedure rational_points(Q)
        // Chabauty
        for p in PrimesInInterval(3, 40) do
            try 
                time data := coleman_data(Q, p, 25);
                time L,v := effective_chabauty(data : bound := 1000, e := 50);
                printf "found %o residue disks.\n", #L;
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
end procedure;

// we know from [ANTS paper] that r < g and differences of rat pts generate a finite index subgroup of J(Q)

// N = 176
Q := y^2 - (x^10 - 8*x^8 - 4*x^7 + 24*x^6 + 16*x^5 - 36*x^4 - 32*x^3 + 16*x^2 + 16*x);
rational_points(Q);

// N = 180
Q := y^2 - (x^2 + x + 1)*(x^4 - 5*x^3 + 12*x^2 - 5*x + 1);
rational_points(Q);

// N = 184
Q := y^2 - (x^6 - 4*x^5 + 8*x^4 - 16*x^3 + 28*x^2 - 24*x + 8);
rational_points(Q);