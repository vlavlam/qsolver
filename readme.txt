Program title: qsolver
Program version: 1.0
Program authors: Vladimir Meshkov, Asen Pashov, Andrey Stolyarov
Corresponding author: Vladimir Meshkov
Contact e-mail: MeshkovVV@my.msu.ru
Developer's repository: https://github.com/vlavlam/qsolver
CPC Program Library archive: https://doi.org/10.17632/x4dgx6k5g8.1
Programming language: Fortran 2003
License: MIT (see the LICENSE file)
The distributed source files qsolver.f90, qpotentials.f90, and qtest.f90
have been validated to conform to the Fortran 2003 standard.

1. Purpose of the package
-------------------------
"qsolver" is a Fortran 2003 program for the highly accurate numerical solution of radial Schrodinger equations for bound and resonant (quasi-bound) states. The program solves both single-channel and coupled-channel (CC) problems, including optionally first-derivative terms, by means of a polynomial collocation method conjugated with the analytical mapping procedure. The "qsolver" supports eigenvalue and eigenfunction calculations of bound states and uses the exterior complex scaling method to search for energies and widths of resonant states. The module optionally provides a straightforward extrapolation of the conventional real-valued potentials and radial coupling functions defined by a user into the complex plane, automatically filters spurious ("ghost") resonance roots, and estimates the first partial derivatives of eigenvalues with respect to model Hamiltonian parameters.

2. Contents of the package
--------------------------
The package contains the files:

  qsolver.f90
      Main source file. Defines module "qsolver_mod" and the public subroutine "qsolver".

  qpotentials.f90
      Module potentials used by the test program (all built-in model potentials and helper functions).

  qtest.f90
      Test driver and worked examples. Compiles into an executable that exercises the main features of "qsolver".

  qtest_out.txt
      Sample output from a comprehensive run of "qtest", demonstrating correct operation.

  LICENSE
      Full text of the MIT License, under which "qsolver" is distributed.

  readme.txt
      This file.

  CITATION.cff
      Citation metadata for the program and its associated publication.

No separate input files are required for the supplied test program: all test potentials and settings are defined directly in qpotentials.f90 and qtest.f90.

3. Documentation
----------------
The complete user manual for the qsolver subroutine is provided as the
extensive header comment block at the top of qsolver.f90. It documents
every keyword argument and the conventions
of the keyword interface, and should be consulted as the primary reference
when calling "qsolver".

The test program qtest.f90 serves two roles: it validates the package
against reference and literature values, and it provides a collection of
worked examples that illustrate how to call "qsolver" for a wide range of
particular problems. When setting up a new calculation, a good starting point is the test closest to the problem at hand.

The associated publication gives the theoretical background and a
discussion of the numerical method used.

4. External requirements
------------------------
To compile and run the package you need:

  - A Fortran 2003 compliant compiler (for example, "gfortran" and "ifort")
  - the external libraries: BLAS, LAPACK and ARPACK

The code has been prepared for standard Unix-like environments (Linux/macOS). On other systems the same source files should work provided equivalent LAPACK/BLAS/ARPACK libraries are available.

Installing the libraries
.........................
On most systems these libraries are available as ready-made packages,
so they do not have to be built from source. Typical commands are:

  Linux (Debian/Ubuntu):
      sudo apt install gfortran liblapack-dev libblas-dev libarpack2-dev

  Linux (Fedora/RHEL):
      sudo dnf install gcc-gfortran lapack-devel blas-devel arpack-devel
      (arpack-devel is provided by the EPEL repository on RHEL/derivatives)

  Linux (Arch):
      sudo pacman -S gcc-fortran lapack blas arpack

  macOS (Homebrew, https://brew.sh):
      brew install gcc lapack openblas arpack
      If the compiler does not find the libraries, point it at the
      Homebrew prefix, for example:
          gfortran ... -L$(brew --prefix)/lib -llapack -lblas -larpack

  Intel compilers (ifort/ifx):
      Intel's oneAPI toolkits ship with the Math Kernel Library (MKL),
      which already provides LAPACK and BLAS. With MKL you only need to
      add ARPACK separately and link MKL instead of -llapack -lblas, e.g.
      on Linux/macOS:
          ifort ... -larpack -qmkl
      (On Windows the corresponding MKL switch is /Qmkl, the Intel
      compilers also require a Microsoft Visual C++ installation, and
      ARPACK needs the extra linking step described in the Windows notes
      below. Both classic ifort and the newer ifx have been used.)

  Windows:
      Two native routes have been verified to build qsolver and pass the
      full qtest suite on Windows (in addition to the Windows Subsystem
      for Linux, WSL, under which the Linux instructions above apply
      directly):

      (a) MSYS2 / MinGW-w64 with gfortran (https://www.msys2.org):
              pacman -S mingw-w64-x86_64-gcc-fortran
                        mingw-w64-x86_64-lapack mingw-w64-x86_64-arpack
          (the arpack package pulls in BLAS/LAPACK). Then use the gfortran
          compilation sequence of Section 5 from a "MINGW64" shell. Keep
          the MinGW-w64 bin directory (e.g. C:\msys64\mingw64\bin) on PATH
          when running qtest, so the BLAS/LAPACK/ARPACK and gfortran
          runtime DLLs are found.

      (b) Intel oneAPI (classic ifort or the newer ifx). The Intel
          compilers additionally require a Microsoft Visual C++
          installation (e.g. the "Build Tools for Visual Studio") for the
          linker and C runtime. Use MKL for BLAS/LAPACK via /Qmkl. A native
          MSVC build of ARPACK is not generally available; the usual ARPACK
          packages (vcpkg arpack-ng, or the MSYS2 package above) are built
          with MinGW gfortran, whose external names are lowercase with a
          trailing underscore. Match that name mangling and link ARPACK's
          import library, for example:
              ifort /names:lowercase /assume:underscore ^
                    qsolver.f90 qpotentials.f90 qtest.f90 /Qmkl ^
                    /link /LIBPATH:"<path-to-arpack>\lib" libarpack.dll.a
          At run time keep the ARPACK DLL and its MinGW dependencies on
          PATH (libarpack, the BLAS/LAPACK DLLs it uses, and the gfortran
          runtime libgfortran, libquadmath, libgcc_s_seh, libwinpthread).

(OpenBLAS is an optimized drop-in replacement for the reference BLAS,
and on some systems also supplies LAPACK; it can be linked in place of
-lblas where available.)

Library home pages and sources:

  BLAS:    https://netlib.org/blas/     (OpenBLAS: https://www.openblas.net)
  LAPACK:  https://netlib.org/lapack/   (reference repository:
           https://github.com/Reference-LAPACK/lapack)
  ARPACK:  the maintained version is ARPACK-NG,
           https://github.com/opencollab/arpack-ng

5. Compilation
--------------
Because qtest.f90 uses modules defined in qsolver.f90 and qpotentials.f90, these files must be compiled first, or all files must be compiled in an order that produces qsolver_mod.mod and potentials.mod before qtest.f90 is compiled.

Verified compilation sequence:

  gfortran -c qsolver.f90
  gfortran -c qpotentials.f90
  gfortran qsolver.o qpotentials.o qtest.f90 -llapack -lblas -larpack -o qtest

If the libraries are not in the default search path, add the appropriate -L option(s), for example:

  gfortran -c qsolver.f90
  gfortran -c qpotentials.f90
  gfortran qsolver.o qpotentials.o qtest.f90 -L/path/to/libs -llapack -lblas -larpack -o qtest

A one-line alternative may also be used if the compiler resolves module dependencies correctly, for example:

  gfortran qsolver.f90 qpotentials.f90 qtest.f90 -llapack -lblas -larpack -o qtest

6. Running the test program
---------------------------
After successful compilation, execute:

  ./qtest

Standard output reports the progress of the tests and a short final summary.
A detailed report is written to:

  qtest_out.txt

The distributed test program contains 32 tests covering:

  - free particle in a box
  - harmonic oscillator
  - eigenvalue problem for the Bessel equation
  - Coulomb potential (hydrogen atom)
  - Morse oscillator
  - double-minimum potentials
  - Lennard-Jones potentials
  - resonant states obtained by complex scaling
  - automatic extrapolation of the user-provided real potentials into the complex plane
  - wavefunctions, expectation values, and energy derivatives
  - derivatives of energies with respect to a mass
  - one- and two-channel benchmark problems with resonant states

7. Expected output from the comprehensive test run
--------------------------------------------------
A successful execution prints a sequence of PASSED messages and ends with a summary similar to:

  Number of passed tests:  32
  Number of failed tests:   0
  Detailed test results are in qtest_out.txt

The detailed file qtest_out.txt contains numerical comparisons against reference or literature values for the included benchmark problems.

8. Using "qsolver" in your own program
------------------------------------
"qsolver" is distributed as a single self-contained module (qsolver.f90)
exposing one public subroutine, rather than as a multi-file library, so it
can be dropped directly into a user's project with no separate library to
build, install, or link. Its only external dependencies are BLAS, LAPACK
and ARPACK. The header comment block of qsolver.f90 documents every
argument and serves as the user manual (see Section 3).

To use "qsolver" in another Fortran program:

  1) compile qsolver.f90
  2) add "use qsolver_mod" in your source
  3) call "qsolver" with keyword arguments

Minimal example (computes the 10 lowest energies of the particle-in-a-box,
or infinite-well, problem on [0, pi]):

  use qsolver_mod
  implicit none
  real(8), parameter :: pi = acos(-1.0d0)
  real(8) :: Ev(10)

  call qsolver(rmin = 0.0d0, &
       &       rmax = pi,    &
       &       mu   = 0.5d0, &
       &       N    = 40,    &
       &       nE   = 10,    &
       &       E    = Ev)

For the full list of arguments and their meanings, see the qsolver.f90
header comment block (Section 3).

9. License
----------
"qsolver" is distributed under the MIT License. The full license text is
provided in the accompanying LICENSE file, and the same license is declared
in the "Licensing provisions" field of the associated publication.

  License: MIT

10. Citation
------------
Associated publication:
  Vladimir Meshkov, Asen Pashov, Andrey Stolyarov,
  "qsolver: A program for the accurate solution of coupled radial
  Schrödinger equations for bound and resonant states",
  Computer Physics Communications 329 (2026), article 110393.
  DOI: https://doi.org/10.1016/j.cpc.2026.110393

If you use qsolver in your research, please cite the publication above.
The article DOI identifies the paper; the CPC Program Library DOI
identifies the archived program files.
The CITATION.cff file provides machine-readable citation metadata, with
the article as the preferred citation for GitHub's "Cite this repository".

BibTeX:

  @article{Meshkov2026qsolver,
    author  = {Meshkov, Vladimir and Pashov, Asen and Stolyarov, Andrey},
    title   = {{qsolver}: A program for the accurate solution of coupled radial
               {Schr{\"o}dinger} equations for bound and resonant states},
    journal = {Computer Physics Communications},
    volume  = {329},
    pages   = {110393},
    year    = {2026},
    doi     = {10.1016/j.cpc.2026.110393}
  }
