# ZeroT-RandomField-JEO
Zero-Temperature Phase Transitions of LiHo<sub>x</sub>Y<sub>1-x</sub>F<sub>4</sub> in the Presence of Random Fields Using Jaded Extremal Optimization

#### Overview
This repository contains the Java source code developed during my undergrad project, titled "Zero-Temperature Phase Transitions of LiHo<sub>x</sub>Y<sub>1-x</sub>F<sub>4</sub> in the Presence of Random Fields." The project explores zero-temprerature phase transitions in disordered magnetic systems. The core focus is on simulating the behavior of the dipolar disordered magnet LiHo<sub>x</sub>Y<sub>1-x</sub>F<sub>4</sub> under varying conditions using a heuristic optimization method called Jaded Extremal Optimization (JEO) \[2,3,7\] to determine the ground states of the system.

#### Repository Contents
- `Main.java`: Main class that orchestrates the simulation runs.
- `ArrayIterator.java`, `Heap.java`: Utility classes that support data management within the simulations.
- `EwaldSum.java`: Implements the Ewald Summation technique for calculating long-range dipolar interactions.
- `PhaseTransCalc.java`: Handles the computation of phase transitions using finite-size scaling.
- `SingleSpin.java`: Manages the properties and behaviors of a single spin within the simulation.
- `trial.java`: Contains trial runs and test cases to verify the correctness of the implemented methods.

#### Running the Simulations
1. Compile the Java files using `javac -cp lib/*;. src/**/*.java -d bin` in the command line.
2. Run the simulation using `java -cp lib/*;bin Main`.

#### Project Background
The simulation investigates the zero-temperature behavior of the LiHo<sub>x</sub>Y<sub>1-x</sub>F<sub>4</sub> compound in a transverse magnetic field, focusing on the phase transition from ferromagnetic to quasi-spin glass states as a function of random field strength. The methodology includes Extremal Optimization for finding ground states and finite-size scaling analysis to determine phase transitions.

#### References
\[1\] Juan Carlos Andresen, Creighton K. Thomas, Helmut G. Katzgraber, and Moshe Schechter. Novel disordering mechanism in ferromagnetic systems with competing interactions. Phys. Rev. Lett., 111: 177202, Oct 2013, https://link.aps.org/doi/10.1103/PhysRevLett.111.177202.

\[2\] S. Boettcher. Extremal optimization for Sherrington-Kirkpatrick spin glasses. The European Physical Journal B - Condensed Matter and Complex Systems, 46(4): 501–505, Aug 2005, https://doi.org/10.1140/epjb/e2005-00280-6.

\[3\] Stefan Boettcher and Allon G. Percus. Optimization with extremal dynamics. Phys. Rev. Lett., 86: 5211–5214, Jun 2001, https://link.aps.org/doi/10.1103/PhysRevLett.86.5211.

\[4\] Michel J P Gingras and Patrik Henelius. Collective phenomena in the LiHo<sub>x</sub>Y<sub>1-x</sub>F<sub>4</sub> quantum ising magnet: Recent progress and open questions. Journal of Physics: Conference Series, 320(1): 012001, 2011, http://stacks.iop.org/1742-6596/320/i=1/a=012001.

\[5\] H.G. Katzgraber. Introduction to Monte Carlo Methods. ArXiv e-prints, May 2009, http://adsabs.harvard.edu/abs/2009arXiv0905.1629K.

\[6\] Helmut G. Katzgraber, Mathias Körner, and A. P. Young. Universality in three-dimensional ising spin glasses: A monte carlo study. Phys. Rev. B, 73: 224432, Jun 2006, https://link.aps.org/doi/10.1103/PhysRevB.73.224432.

\[7\] A. Alan Middleton. Improved extremal optimization for the Ising spin glass. Phys. Rev. E, 69: 055701, May 2004, https://link.aps.org/doi/10.1103/PhysRevE.69.055701.

\[8\] Moshe Schechter. LiHo<sub>x</sub>Y<sub>1-x</sub>F<sub>4</sub>. Phys. Rev. B, 77: 020401, Jan 2008, https://link.aps.org/doi/10.1103/PhysRevB.77.020401.

\[9\] P. Stasiak and M. J. P. Gingras. Evidence for a Finite-Temperature Spin Glass Transition in a Diluted Dipolar Heisenberg Model in Three Dimensions, Appendix C. ArXiv e-prints, December 2009, https://arxiv.org/abs/0912.3469.

\[10\] Moshe Schechter and Nicolas Laflorencie. Quantum spin glass and the dipolar interaction. Phys. Rev. Lett., 97: 137204, Sep 2006, https://link.aps.org/doi/10.1103/PhysRevLett.97.137204.
