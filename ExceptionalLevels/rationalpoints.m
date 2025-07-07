load "../Coleman/coleman.m";
//40, 48, 72, 80, 88, 96, 99, 100, 108, 112, 120, 135, 144, 147, 162, 176, 180, 184, 196, 200, 216, 224, 225, 240
N := 40; //2^3*5
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//Elliptic Curve defined by y^2 = x^3 + x^2 + 4*x + 4 over Rational Field
//Canonical embedding: yes
SIntegralPoints(X, [2]);
// [ (-1 : 0 : 1), (0 : -2 : 1), (4 : -10 : 1) ]
//only one given up to hyperelliptic involution and no infinite point
SIntegralPoints(X, []);
//[ (-1 : 0 : 1), (0 : -2 : 1), (4 : -10 : 1) ]
Rank(X);
// 0 true
TorsionSubgroup(X);
// Abelian Group isomorphic to Z/6
// Defined on 1 generator
// Relations:
// 6*$.1 = 0
//6 rational points


N := 48; //2^3*3
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//Elliptic Curve defined by y^2 = x^3 - x^2 - 4*x + 4 over Rational Field
//Canonical embedding: yes
SIntegralPoints(X, [2]);
// [ (-2 : 0 : 1), (0 : -2 : 1), (1 : 0 : 1), (2 : 0 : 1), (4 : 6 : 1) ]
SIntegralPoints(X, []);
//[ (-2 : 0 : 1), (0 : -2 : 1), (1 : 0 : 1), (2 : 0 : 1), (4 : 6 : 1) ]
Rank(X);
// 0 true
TorsionSubgroup(X);
// Abelian Group isomorphic to Z/2 + Z/4
// Defined on 2 generators
// Relations:
// 2*$.1 = 0
// 4*$.2 = 0
// 8 rational points


N := 72; //2^3*3^2 
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
//canonical embedding: yes
SIntegralPoints(X, [2,3]);
// [ (-1 : 0 : 1), (0 : -1 : 1), (2 : -3 : 1) ]
Rank(X);
//0 true
TorsionSubgroup(X);
// Abelian Group isomorphic to Z/6
// Defined on 1 generator
// Relations:
// 6*$.1 = 0
// 6 rational points


N := 80; //2^4*5
//Magma gives wrong model here
X := EllipticCurve([0, 1, 0, -1, 0]);
//Elliptic Curve defined by y^2 = x^3 + x^2 - x  over Rational Field;
//canonical embedding: yes
SIntegralPoints(X, [2]);
// [ (-1 : 0 : 1), (0 : -2 : 1), (4 : -10 : 1) ]
Rank(X);
//0 true
TorsionSubgroup(X);
// Abelian Group isomorphic to Z/6
// Defined on 1 generator
// Relations:
// 6*$.1 = 0s
// 6 rational points

N := 88;
f := x^6 + 2*x^5 - x^4 + 12*x^3 - x^2 + 2*x + 1;
H := HyperellipticCurve(f);
rl, ru, gens := RankBounds(f,2 :ReturnGenerators := true);
assert rl eq ru and rl eq 1;
chab, s1, s2 := Chabauty(gens[1]);
chab;
//{ (0 : -1 : 1), (1 : -4 : 1), (1 : 4 : 1), (0 : 1 : 1), (1 : -1 : 0), (1 : 1 : 0) }
// 6 rational points

N := 96; //2^4 *3
//Magma gives wrong model here
X := EllipticCurve([0, -1, 0, 1, 0]);
//Elliptic Curve defined by y^2 =  x^3 - x^2 + x over Rational Field
//cannonical embedding: yes
SIntegralPoints(X, [2]);
// [ (-2 : 0 : 1), (0 : -2 : 1), (1 : 0 : 1), (2 : 0 : 1), (4 : 6 : 1) ]
Rank(X);
// 0 true
TorsionSubgroup(X);
// Abelian Group isomorphic to Z/4
// Defined on 1 generator
// Relations:
// 4*$.1 = 0
// 4 Rational points

N := 99; //3^2*11
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 - 2*x over Rational Field
//cannonical embedding: yes
SIntegralPoints(X, [3]);
// [ (-1 : 0 : 1), (-8/9 : 13/27 : 1), (0 : 0 : 1), (2 : -3 : 1), (3 : 2 : 1), (26 : -144 : 1) ]
Rank(X);
//1 true
// INFINITE RATIONAL POINTS


N := 100; //2^2*5^2
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 3*x + 1 over Rational Field
//cannonical embedding: yes
SIntegralPoints(X, [2,5]);
// [ (-1 : 2 : 1), (1 : 0 : 1) ]
Rank(X);
// 0 true
TorsionSubgroup(X);
// Abelian Group isomorphic to Z/5
// Defined on 1 generator
// Relations:
// 5*$.1 = 0
// 5 rational points


N := 108;
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//Elliptic Curve defined by y^2 + x*y + y = x^3 - x^2 + x - 1 over Rational Field
//cannonical embedding: yes
SIntegralPoints(X, [2,3]);
// [ (1 : 0 : 1) ]
Rank(X);
// 0 true
TorsionSubgroup(X);
// Abelian Group isomorphic to Z/3
// Defined on 1 generator
// Relations:
// 3*$.1 = 0



N := 112;
f := x^6 - 2*x^5 + 11*x^4 - 4*x^3 + 11*x^2 - 2*x + 1 ;
rl, ru, gens := RankBounds(f,2 :ReturnGenerators := true);
assert rl eq ru and rl eq 1;
chab, s1, s2 := Chabauty(gens[1]);
//[(pi^(-1))(P) : P in chab];
chab;
//[ (1 : 1 : 1), (0 : 2 : 1), (1 : 0 : 0), (1 : 1 : 0), (1 : 0 : 1), (0 : -2 : 1) ]
// 6 rational points


N := 120;
X :=  EllipticCurve([0, 1, 0, -1, 0]);
//Model is wrong in magma
//Elliptic Curve defined by y^2 =  x^3 + x^2 - x  over Rational Field
//canonical emb: yes
SIntegralPoints(X, [2,3]);
// [ (-1 : 0 : 1), (0 : -2 : 1), (4 : -10 : 1) ]
Rank(X);
// 0 true
// TorsionSubgroup(X);
// Abelian Group isomorphic to Z/6
// Defined on 1 generator
// Relations:
// 6*$.1 = 0
//6 rational points


N := 135;
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//genus 2
C, pi := SimplifiedModel(X);
//Hyperelliptic Curve defined by y^2 = x^6 - 4*x^5 + 10*x^4 - 20*x^3 + 17*x^2 - 20*x over Rational Field
f, h := HyperellipticPolynomials(C);
rl, ru, gens := RankBounds(f,2 :ReturnGenerators := true);
assert rl eq ru and rl eq 1;
chab, s1, s2 := Chabauty(gens[1]);
[(pi^(-1))(P) : P in chab];
//[ (0 : 0 : 1), (1 : 0 : 0), (1 : 1 : 0) ]
//3 rational points


N := 144; //this is rank 0, so we can just run Chabauty
//-y^3 + (-x^2 + 4*x - 5)*y^2 + (2*x^3 - 6*x^2 + 4*x - 1)*y - x^4 + 2*x^3 - x^2
//genus 3
Q:= y^3 - y^2*x^2 - 4*y^2 + 2*y*x^3 - 4*y*x^2 + 4*y*x + 6*y - x^4 + 4*x^3- 8*x^2 + 8*x - 7;
data := coleman_data(Q, 13, 20);
L,v :=effective_chabauty(data:bound:=1000,e:=50); 
rpts := Q_points(data,1000);
#L eq #rpts;
// {@ (0 : 1 : 0), (0 : 0 : 1), (1 : 1 : 0), (1 : 0 : 1), (1 : -1 : 1), (-1 : -1 : 1) @}
//6 rational points

N := 147;
//already done in https://arxiv.org/abs/2203.05541, p.16
//Hyperelliptic Curve defined by y^2 = x^6 + 6*x^5 + 11*x^4 + 6*x^3 + 5*x^2 + 4 
//genus 2
//12 rational points


N := 162;
Q := y^3 - y^2*x^2 - 3*y^2 + 3*y*x^2 - 3*y*x + 3*y + 2*x^3 - 3*x^2 + 3*x;
//genus 3
rpts := Q_points(data, 1000);
data := coleman_data(Q, 5, 10);
Qpoints := Q_points(data,1000);
P1 := Qpoints[1];
P2 := Qpoints[2];
IP1P2, N2 := coleman_integrals_on_basis(P1, P2, data: e:=50); //verify that we have a non-torsion rational point
assert not IsWeaklyZero(IP1P2[1]);
L,v :=effective_chabauty(data:bound:=1000,e:=50); 
#L eq #rpts;
// {@ (-1/2 : 1/2 : 1), (0 : 1 : 0), (0 : 0 : 1), (1 : 0 : 0), (1 : 1/2 : 1) @}
//5 rational points


N := 196; //this is rank 0, so we can just run Chabauty
Q := y^3 - y^2*x^2 - 2*y^2 + 4*y*x^2 - y*x* - 2*y - x^3;
data := coleman_data(Q, 5, 10);
L,v :=effective_chabauty(data:bound:=1000,e:=50); 
rpts := RationalPoints(X :Bound := 1000);
#L eq #rpts;
//{@ (0 : -1 : 1), (0 : 0 : 1), (1 : 1 : 0), (-1 : 1 : 0) @}
//4 rational points


N := 240; //This is rank 1
X := X0NQuotient(N, [p^Valuation(N,p) : p in PrimeDivisors(N)]);
//genus 3
Q:= y^3 - y^2*x^2 - y^2*1^2 + 2*y*x^3 + 2*y*x^2+ y - x^4 - 2*x^3 - x^2;
data := coleman_data(Q, 11, 10);
Qpoints := Q_points(data,1000);
P1 := Qpoints[1];
P2 := Qpoints[2];
IP1P2, N2 := coleman_integrals_on_basis(P1, P2, data: e:=50); //verify that we have a non-torsion rational point
L,v :=effective_chabauty(data:bound:=1000,e:=50); 
#L eq #Qpoints;
// {@ (1 : 1 : 1), (0 : 1 : 0), (0 : 0 : 1), (-1 : 0 : 1), (-1 : 1 : 1), (1 : 1 :0) @}
//6 rational points



