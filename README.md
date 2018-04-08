# GSCheckDuplicatesCharge

The comparison between two structures is based on the proposal by  Grigoryan et al.(1) and Rogan et al.(2). The criteria for
uniqueness of the geometry are based on mulliken atomic charges (g03 and g09 outputs files) (3).

Getting started

  a) - Prerequisites:

    GSCheckDuplicates is written in perl. The program has only been tested on Mac OS X, Linux and Windows.

  b) - Running GSCheckDuplicates program:

    perl GS_Similarity.pl 0.05 

    * 0.05 is threshold duplicate.

    * The program call data base ./Outs/






(1) Grigoryan, V. G.; Springborg, M. Structure and Energetics of Ni Clusters with up to 150 Atoms. Chem. Phys. Lett. 2003, 375 (1–2), 219–226.

(2) Rogan, J.; Ramirez, M.; Varas, A.; Kiwi, M. How Relevant Is the Choice of Classical Potentials in Finding Minimal Energy Cluster Conformations? Comput. Theor. Chem. 2013, 1021, 155–163.

(3) http://gaussian.com/g09citation/
