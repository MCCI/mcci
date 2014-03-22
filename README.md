Monte Carlo Configuration Interaction
=====================================

Welcome to the web page of the Monte Carlo Configuration Interaction (MCCI)
program. MCCI calculates electronic structure of molecular systems using a
method based on Configuration Interaction (CI). In conventional CI, all
configurations (usually spin-projected Slater Determinants, or Configuration
State Functions (CSFs)) up to a given excitation level are included in the wave
function. While this produces a highly accurate wave function, the number of
CSFs to include (and thus also the dimensions of the matrix to be
diagonalized) increases combinatorially with the truncation level. In addition,
the fixed truncation level of conventional CI causes the method to be
non-size-consistent.
MCCI is based on the observation that in a given vector from a truncated CI
run, only a small fraction of the CSFs produce a significant contribution to
the wave function, i.e. the absolute value of their coefficient in the
expansion is small. It follows that the wave functions (and their energies)
can be described to a satisfactory level of accuracy by omitting these
insignificant CSFs from the expansion, which keeps the problem computationally
tractable. In MCCI, the truncation of the expansion is determined solely by
the contributions to the wave function of the individual CSFs; in an iterative
procedure, new CSFs are generated from the CSFs in the current vector, and
either kept or discarded depending on their coefficient in the wave function.
Note that in this fashion, arbitrarily highly excited CSFs will be included (as
long as their coefficient remains above the threshold), so in principle, size
consistency is achieved.
For more background information on the formalism behind MCCI, please see the
first two papers in the section "Publications".

As its input data, MCCI requires a list of molecular orbital (MO) operator
integrals of the Hamiltonian's one- and two-body operators. Any electronic
structure program that can output these integrals can in principle be used.

Examples
========
Since the mcci algorithm relies on random sampling of the CSF space and retains
CSFs based on their contribution to the wave function, it encounters no
problems describing systems with significant multireference (i.e. more than one
significant CSF) character. Even with a large significance threshold (e.g.
cmin=5.0 E-3), the algorithm will reliably find the most significant CSFs, at a
significantly reduced computational cost compared to other multi-configuration
methods.
Also of note is the scaling behaviour of the code; the figure below shows a
graphical interpretation of the CI matrix obtained from a run of mcci. The
stronger interactions (in red) all involve a limited set of CSFs, which we call
the reference set. This implies that once this reference space is captured in
the exapansion, different nodes can sample the extended space independently
(since those CSFs won't interact strongly with each other), producing
near-linear scaling.

Compile Instructions
====================

The provided Makefile assumes the availability of an MPI-enabled
Fortran 90 compiler under the name 'mpif90'. This can be changed
manually, but the user needs to check that all the flags used in
the Makefile are supported by the other compiler.

If the compiler is changed, please note that the system clock
called in the subroutine gen_seed.f90 may also need to be
changed, depending on the behaviour of the system_clock standard
subroutine on the platform used.


Keywords for the mcci.in file
-----------------------------

    ieig              =  Eigenvalue to be calculated. Each eigenvalue in each irrep begins with ieig=1 
                         i.e. ieig=1 is lowest energy eigenvalue in an irrep, ieig=2 is first excitation in an irrep, ... 
    n_up              =  Number of spin up (alpha) electrons
    n_dn              =  Number of spin down (beta) electrons
    cmin              =  Coefficient threshold:  Defines the minimum value of CSF's coefficient to be kept following
                         a matrix diagonalization, i.e. the "pruning threshold"
    s                 =  Spin ( units of h_bar ) (principle case M_s=S)
    maxtry            =  The maximum number of cycles (branch, diagonalization, prune) that should be performed
    mo_up(:)          =  Molecular orbital indices for up (alpha) electrons
    mo_dn(:)          =  Molecular orbital indices for down (beta) electrons
    lmin              =  Minimum vector length- branch will always generate lmin CSFs
    npfull            =  A "full prune" step will be performed every npfull cycles
    lkeep             =  Configurations to keep throughout a calculation. The first lkeep CSFs will not be pruned.
    hmin              =  H matrix threshold- elements of the H matrix with absolute values below hmin are not stored
    davidson_stop     =  Maximum number of iterations before stopping for convergence in the Davidson algorithm
    bmin              =  Minimum vector boost (used to adjust Monte Carlo sample size in early stages of a calculation)
    bmax              =  Maximum vector boost (used to adjust Monte Carlo sample size in early stages of a calculation)
    frac              =  Branching ratio: number of new CSFs / number of CSFs in current CI vector
    conv_thresh_e     =  Energy convergence: change in energy between npfull cycles is less that conv_thresh_e after ncheck tests
    conv_thresh_l     =  Vector convergencer: change in vector length between npfull cycles is less than conv_thresh_l for ncheck tests 
                         (Note, the number of tests, n, is defined by the keyword conv_history)
    restart           = .TRUE. or .FALSE. If this is set to true, then, a calculation will be 
                         restarted from a previous calculation and civ_in must be present (i.e. previous runs civ_out file 
                         should be renamed to civ_in).  If set to false, a calculation will begin from the CSF with
                         occupations define in mcci.in
    test              = .TRUE. or .FALSE. Debugging flag
    time              = .TRUE. or .FALSE. Timing information, general
    time_all          = .TRUE. or .FALSE. Timing information, detailed
    nobrnch_first     = .TRUE. or .FALSE. No branching in the first step
    nodiag            = .TRUE. or .FALSE. Used for collecting CSFs only (No matrix diagonalisation is peformed during cycles)
    i_want_conv       = .TRUE. or .FALSE. Stop mcci by using convergence criteria.
    npfull_conv       = .TRUE. or .FALSE. Convergence test only in npfull steps
    conv_average      =  Number of npfull steps to be averaged for convergence
    conv_history      =  Window of npfull steps to be include in convergence checking
    frozen_doubly     =  Indices of orbitals to be frozen (i.e. no excitations are allowed from these molecular orbitals 
                         during the CI calculation). Note frozen orbitals must be doubly occupied.
    mo_active         =  Orbital lists of active orbitals, includes all occupied and (initially) unoccupied orbitals to include 
                         in the CI calculations. Note excludes frozen_doubly and all inactive virtual orbitals
    lref              =  Branch (i.e. generate new CSFs) relative to only the first lref CSFs in a CI vector


Example: setting up a calculation in a higher symmetry
======================================================

From the original control file that is obtained from the geometry relaxed
system, note the numbers for the occupied symmetries, e.g.

    A1 1-7
    A2 1-3
    B1 1-9
    B2 1-2

Then run eiger -a and go as high up in energy as required and note the numbers
of orbitals that you want to run for mpgrad. Make sure to leave no gaps
between representations.

Then change A1 - B2 for the new symmetries and run mpgrad_ints twice.

To fill this into mcci, we only need to look at the latest version of the
control file.

Fill in the n_up and n_down as required. Then we must label the symmetry
representations in correct order, so if have:

    MPGRAD control     MCCI counting
    
    A1 1-10            1-10
    A2 1-5             11-16
    B1 1-15            17-...
    B2 1-4             .......

The general rule is: for the first representation it stays as is, for the
following :

(#orbitals) + (#electrons) -1

And when filling in the mcci.in file , the same rule applies.

Publications
============

Method
------

* Estimating Full Configuration Interaction Limits from a Monte Carlo Selection of the Expansion Space, J.C. Greer, Journal of Chemical Physics, 103 pp. 1821-1828 (1995)

* Monte Carlo Configuration Interaction, J. C. Greer, Journal of Computational Physics, 146 pp. 181-202 (1998)

* A Parallelization Model for Successive Approximations to the Rayleigh-Ritz Linear Variational Problem, J.C. Greer, IEEE Transactions on Parallel and Distributed Systems, 9 pp. 938-951 (1998)

* A Monte Carlo Configuration Generation Computer Program for the Calculation of Electronic States of Atoms, Molecules and Quantum Dots, L. Tong, M. Nolan, T. Cheng and J. C. Greer, Computer Physics Communications, 131 pp. 142-163 (2000)

Applications
------------

* Consistent Treatment of Correlation Effects in Molecular Dissociation Studies Using Randomly Chosen Configurations, J. C. Greer, Journal of Chemical Physics, 103 pp. 7996-8003 (1995)

* Electronic Correlation Energy in Linear and Cyclic Carbon Tetramers, J. C. Greer, Chemical Physics Letters, 306 pp. 197-201 (1999)

* Monte Carlo Configuration Interaction Predictions for the Electronic Spectra of Ne, CH2, C2, N2, and H2O Compared to Full Configuration Interaction Calculations, W. Gyoffry, R. J. Bartlett and J. C. Greer, Journal of Chemical Physics 129, 064103 (2008)

* Spin-polarisation Mechanisms of the Nitrogen-Vacancy Center in Diamond, P. Delaney, J.C. Greer, and J.A. Larsson, Nano Letters, 10, 610-614 (2010)

Miscellaneous
=============

*  Impact of Electron-Electron Cusp on Configuration Interaction Energies, D. Prendergast, M. Nolan, C. Filippi, S. Fahy, and J. C. Greer, Journal of Chemical Physics, 115 pp. 1626-1634 (2001)

* Tools for Analysing Configuration Interaction Wavefunctions, P. Delaney and J. C. Greer, Computational Materials Science, 28 pp. 240-249 (2003)

* Statistical Estimates of Electron Correlations, W. Gyoffry, T. M. Henderson and J. C. Greer, Journal of Chemical Physics, 125, 054104 (2006)

* Determining Complex Absorbing Potentials from Electron Self Energies, T. M. Henderson, G. Fagas, E. Hyde, and J. C. Greer, Journal of Chemical Physics, 125,  244104 (2006)

Applications in Molecular Electronics (with the program VICI)
=============================================================

* Correlated Electron Transport in Molecular Electronics, P. Delaney and J. C. Greer, Physical Review Letters, 93 pp. 036805-036808 (2004)

* Quantum Electronic Transport in a Configuration Interaction Basis, P. Delaney and J. C. Greer, International Journal of Quantum Chemistry, 100 pp.1163-1169 (2004)

* Independent Particle Descriptions of Tunneling Using the Many-body Quantum Transport Approach, G. Fagas, P. Delaney and J. C. Greer, Physical Review B, 73, 241314(R) (2006).

* Tunnelling in Alkanes Anchored to Gold Electrodes via Amine Groups, G. Fagas and J. C. Greer, Nanotechnology, 18 424010 (2007)

* A Comparative Study for the Calculation of Tunnel Currents Across Silane Diamines/Dithiols and Alkane Diamines/Dithiols, Shane McDermott, Chris George, Giorgos Fagas, J. C. Greer, and M. A. Ratner, Journal of Physical Chemistry C, 113, pp. 744-750 (2009) 

Contributors
============
* Jim Greer
* Werner Gyorffy
* Thomas Kelly
* Mark Szepieniec
* Jeremy Coe

Contact
=======

For questions, comments or troubleshooting assistance, please email jim.greer@tyndall.ie
