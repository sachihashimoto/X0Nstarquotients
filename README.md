# X0Nstarquotients
Code associated to the paper "Rational points on $X_0(N)^*$ when $N$ is non-squarefree" by Sachi Hashimoto, Timo Keller, and Samuel Le Fourn

The code requires Magma and Sage, and has been tested on Magma V2.28-18, Magma V2.28-13 and Sage 10.5, 10.6.

The files in the main folder are as follows:

- gonal_maps.m is modified Magma code to output several gonal maps instead of only one

- heegner.m computes Heegner points and their associated data (such as admissible prime ideals) on the star quotient. The main functions are HeegnerPoints and galoisALcompatibleHps

- J0wplusminus.m computes subvarieties of the Jacobian J0(N) with prescribed Atkin-Lehner eigenvalues and their associated rank

- modelsX0Nstar.m is most useful for constructing models of hyperelliptic X0(N)^* 

- Qcusps.m is a short script to compute the number of cusps that are rational on X0(N)^*

- rootsofunity.m computes sums of roots of unity and their support, as described in section 5 of the paper.

- verify.sage is a sage script to search the LMFDB for newforms of certain small levels and specified Atkin--Lehner signs, to check the existence of rank zero quotients of J0(pq) (with prescribed signs)

- X0NstarCCsolver.m is a script to automate the computation of rational points on X0(N)^* by first computing a plane model that can be used in the Chabauty--Coleman code of Balakrishnan and Tuitman, then searching for rational points  checking that the differences known small rational points generate a finite index subgroup of the Mordell--Weil group, and trying small good primes until it (hopefully) succeeds

- X0Nstarmodelgenus1.m produces models of genus 1 star quotients 


The folders are as follows:

- ExceptionalLevels contains all the code and log files associated to computing the list of exceptional levels and classifying their rational points

- SmallLevelComputations contains all the code and log files associated to computing the list of X0(N)^* of genus between 1 and 5 and classifying their rational points

- Coleman is a copy of the repository of Jennifer Balakrishnan and Jan Tuitman for computing Coleman integrals https://github.com/jtuitman/Coleman

- QuadraticPoints is a copy of the repository of Nikola Adzaga, Timo Keller, Philippe Michaud-Jacobs, Filip Najman, Ekin Ozman, and Borna Vukorepa https://github.com/timokellermath/quadraticpoints
