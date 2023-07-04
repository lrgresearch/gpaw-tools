---
layout: default
nav_order: 3
title: Release Notes
parent: development
---

# Release notes

### Development version
* n/a

### Version 23.7.0

* Using dtype as default for PW calculations.
* Drawn DOS and Band figures are aligned with respect to Fermi energy level.
* Implementing autoscale in y-direction for DOS graphs.
* Correcting all possible errors due to gpts and kpts density usage.
* Lowering memory consumption in optical calculations.
* Some examples are changed (Tetrahedron method is used in CrO2 calculations, Bulk-Al calculations are now Bulk-GaAs).
* Using world.size for nblocks in optical RPA calc.
* Updated some calculation default values in both `gpawsolve.py` and `gg.py`.
* Small corrections are done to EXX-related parts, XYYY formatted band result data output, optical RPA calculation, ...etc...
* Using `DielectricFunction` with frequencies.
* PW-EXX mode is removed. EXX can be used directly under PW.
* Mostly, no need to use `outdir` variable. However, still can be used.
* Phonon calculation feature is added to `gg.py`.
* Phonopy version information can be viewed when using -v argument.
* Energy consumption measurement with -e argument. This feature only works with Intel CPUs after Sandy Bridge generation. Results are given in kWh!
* Restart -r, --restart argument is now depreceted. There is nothing to be restarted. Instead of -r argument, the new keyword `Ground_calc` is introduced.
* Basic phonon dispersion calculation feature with Phonopy! At least it works for bulk Al. It is not very mature so please keep this in mind. (Thanks to Michael Lamparski for his help and MIT licensed code that he shared.)
* Save figures in higher dpi.
* Fix some bugs, add new variable and rearrange some variable in `gg.py`.
* Add Energy_min variable. Energy_min and Energy_max variables are now working on both band and DOS figures.
* Remove unnecessary variables in the input file for the examples/Si-2atoms-optical example.
* Band structure data can be now exported in XYYYY type ASCII file (Thanks to Andrej Kesely). 
* Fixed some unused imports and local variables.
* Small fixes at `gg.py` file.
* Adding `struct_from_file()` function for future usage.
* `gpawsolve.py` is not only a Python script anymore. The structure of the `gpawsolve.py` is rewritten. The calculations are related to a class named `gpawsolve`. Also, there are functions related to all possible calculations in this class as: `structurecalc()`, `groundcalc()`, `elasticcalc()`, `doscalc()`, `bandcalc()`, `densitycalc()` and `opticalcalc()`. The structure is still primitive and code is rely mostly on global variables, however it is a start and it will be easy to use when it finished properly.
* From this release "-o" argument is deprecated. Code, all examples and related BASH script are fixed.

### Version 23.2.0

* DOS calculation part is changed completely. All calculations for DOS, PDOS and RawPDOS are done with RawPDOS.
* New variable is added: `DOS_convergence`.
* The variables used in `asapsolve.py` is also changed. From this version, `asapsolve.py` will use a special Snake_case (first letter of each variable is capitalized) variables for the input file variables usage.Variable that are affected are: *Manualpbc -> Manual_PBC, pbcmanual -> PBC_constraints, PotentialUsed -> OpenKIM_potential, SolveDoubleElementProblem -> Solve_double_element_problem.*
* Correction of the creating output file that writes spacegroup and special points to a wrong folder.
* Making `do_all_examples.sh` script executable.
* Wrong usage of GW calculation type in `gg.py` is corrected.
* All new variable changes are corrected and new variables are added to `gg.py`.
* `simple_benchmark2021.py` is simplified and renamed as `simple_benchmark2023.py`.
* With time, `gpawsolve.py` uses many variables in it and in input files. However, it does not have a proper naming convention for these variables. Also some of the variable names are misleading the user. From this version, `gpawsolve.py` will use a special Snake_case (first letter of each variable is capitalized) variables for the input file variables usage. Variables that are affected: *mode -> Mode, fmaxval -> Max_F_tolerance, whichstrain -> Relax_cell, cut_off_energy -> Cut_off_energy, kpts_density -> Ground_kpts_dens, kpts_x -> Ground_kpts_x, kpts_y -> Ground_kpts_y, kpts_z -> Ground_kpts_y, gpts_density -> Ground_gpts_dens, gpts_x -> Ground_gpts_x, gpts_y -> Ground_gpts_y, gpts_z -> Ground_gpts_z, Hubbard -> Setup_params, band_path -> Band_path, band_npoints -> Band_npoints, gridref -> Refine_grid energy_max -> Energy_max, GWtype -> GW_type, GWkpoints -> GW_kpoints_list, GWtruncation -> GW_truncation, GWcut_off_energy -> GW_cut_off_energy, GWbandVB -> GW_valence_band_no, GWbandCB -> GW_conduction_band_no, GWppa -> GW_PPA, GWq0correction -> GW_q0_correction, GWnblock -> GW_nblocks_max, GWbandinterpolation -> GW_interpolate_band, opttype -> Opt_calc_type, optshift -> Opt_shift_en, optBSEvb -> Opt_BSE_valence, optBSEcb -> Opt_BSE_conduction, optBSEminEn -> Opt_BSE_min_en, optBSEmaxEn -> Opt_BSE_max_en, optBSEnumdata -> Opt_BSE_num_of_data, num_of_bands -> Opt_num_of_bands, optFDsmear -> Opt_FD_smearing, opteta -> Opt_eta, optdomega0 -> Opt_domega0, optomega2 -> Opt_omega2, optecut -> Opt_cut_of_energy, optnblocks -> Opt_nblocks, MPIcores -> MPI_cores*
* `find3Dmin.py` A script which draws 3D contour plot of E vs. latticeparams and show the minimum datapoint using the optimize_latticeparam.py's output, is added.

### Version 22.7.0

* `optimize_latticeparam.py` now can work for both lattice params a and c. Also draws 3D fig of Energy dependent latt_a - latt_c.
*  `quickoptimize.py` works like `gpawsolve.py` now. Its name is now `asapsolve.py`. 
* New default optimizer is QuasiNewton (BFGSLineSearch).
* New keyword `Optimizer`. Users can now choose QuasiNewton, GPMin, LFBGS or FIRE minimizer for geometry optimization. 
* `gg.py` includes all new keywords.
* Grid point density or manual grid points for axis (LCAO only).
* Include new keywords for LBFGS geometry optimization `Damping`, `Alpha` and `Max_step`.
* `Geo_optim` keyword for better optimization usage with filters (whichstrain).
* Examples are simplified. Most of the unnecessary keywords are removed.
* Proper logging for LCAO ground state calculations
* Fix LCAO spinpol calculation (Thanks to Toma Susi).
* Include new keyword `Mixer_type`.
* Fix help description text width problem.
* Execution timing data of all calculations are saved to `FILENAME-6-Result-Log-Timings.txt` file.
* Instead of direct execution, all tasks are added to task-spooler queue in `do_all_examples.sh`script.

### Version 22.5.0

* Successful HSE03, HSE06 calculations for ground state, DOS and band structure.
* New example for HSE06 calculations: `Si-with-HSE`.
* Colorize errors, warnings and information output with ANSI codes.
* Proper error handling for restart mode "file not found" error. 
* New keyword `Fix_symmetry` added to `gpawsolve.py`, `gg.py` for preserving the spacegroup symmetry during optimisation.
* Small changes in the `gg.py`
* No need to import ASE object inside optimization scripts. Optimizations are working with CIF files only.
* `gpawsolve.py` now prints previous and final spacegroup information and usable special points information for band structure calculations.
* `-v` argument now shows version information of gpaw-tools, and used GPAW, and ASE. It gives more choice for possible tarball and zipball packages. Also, it does not give an error in case of no internet connection available.
* 3 new keywords `Ground_convergence`, `Band_convergence` and `Occupation` are added to `gpawsolve.py`, `gg.py` and examples.
* Fix `do_all_examples.sh` Bash script.
* Optimization scripts do not need ASE object insertation. They can run with using CIF file as an argument.
* RawPDOS, which gives PDOS over orbitals, is added.
* For band calculations, result file in JSON format is added. This file can be opened with `ase band-structure` command.

### Version 22.4.0

* New optimization script `optimize_kptsdensity.py` for k-point density optimization instead of k-point number optimization.
* `optimize_cutoff.py`, `optimize_kpoints.py` and `optimize_latticeparam.py` have `xc_used` in parameters list.
* Naming of some of the output files are fixed. 
* Bethe-Salpeter Equation (BSE) solution is added to optical calculations.
* 7 new keywords are added for BSE calculations. `opttype`, `optshift`, `optBSEvb`, `optbsecb`, `optBSEminEn`, `optBSEmaxEn`, `optbsenumdata`
* `Si-2atoms-optical` example is now running for both RPA and BSE. Previously, its calculation has 2 steps , now 3 steps.
* CONTRIBUTING and CODE_OF_CONDUCT information is added.
* Fix: Show atom numbers starting from 1 not 0.

### Version 22.3.0

* New calculation: with `Elastic_calc=True`, Equation of State and elastic tensor values will be calculated.
* `Elastic_calc` is added to `gg.py`.
* A new example about elastic calculations is added.
* Many comments added to `gg.py` for better understanding the code.
* `DOS_npoints` and `DOS_width` variables are added for number of points and smearing value, respectively.
* Saving PNG versions of band structure and DOS even in -d argument is not used.
* `shrinkgpw.py`command for extracting wavefunctions from huge gpw files and save with a different name.
* New benchmarks were added.

### Version 21.12.0

* EXX mode is renamed as PW-EXX.
* Default values of variables are changed.
* Previous -i (input) argument is changed as -g (geometry). It is more logical, because it is used for geometry.
* Previous -c (config) argument is changed as -i (input). It is more general, convenient and understandable.
* `gg.py` can now work with just enough number of variables. Previously, it must see all variables.
* `kpts_density` is added.
* Units used in cube file are changed for Bader analysis.
* There is no general config file anymore.
* Small bugfixes.

### Version 21.11.0

* PDOS calculations.
* PW mode can use GLLB-SC xc now.
* `optimize_cutoff.py` and `optimize_kpoints.py` can use CIF, XYZ, etc... files as input file. No need to include the ASE object inside these scripts anymore.
* The nbands parameter is changed in PW mode.
* In GW calculations, calculation could not be done because of interpolation in drawing the figure when the data did not have a minimum of 3 points. Now there is a variable to use interpolation or not.
* In GW calculations, `gpawsolve.py` can write quasiparticle energies to a file separately.
* DFT+U calculation ability is added for PW and LCAO modes.
* `gg.py` is better now. It is compatible with new arguments and removed variables. It can run in any directory and handling of `GWkpoints` and `GWtruncation` variables are correct now.
* A new argument '-d' is added. This argument makes script draw the DOS and band calculation results. In the past, it was a varible in the config file.
* `WantedCIFexport` variable is removed. 
* A new argument '-r' is added. This argument makes script pass the ground state calculations and continue with the next calculation.

### Version 21.10.1

* An important bug made it impossible to work with existing examples with `gg.py`. It is now resolved.

### Version 21.10.0

* 5 different examples are added to show simply different usage cases.
* Initializing magnetic moment problem is solved.
* Version argument is added.
* GW parameters are also added to `gg.py`
* Add some optical parameters to config files and `gg.py`.
* Major change: `gpawsolve.py` and `gg.py` are now working as commands. Config files can be used as general input files, you can put ASE Atom object plus every parameter that `gpawsolve.py` accepts. If you want, you can provide a CIF file instead of using an atom object. You can run `gpawsolve.py` and `gg.py` from any folder.
* Optical: Refractive index, extinction index, absorption, and reflectivity calculations.
* `gg.py` is now opening ase gui when the user clicks the structure image.
* New argument parsing scheme for better future usages.
* Very basic PW-EXX mode with HSE06 and PBE06. (Only some ground-state calculations.)
* Adding GW0 and G0W0-GW0 selector.
* Adding GW approximation to `gpawsolve.py` (only bands).
* Many other small corrections.

### Version 21.9.0

* Corrected `quickoptimize.py` behaviour.
* Many code quality and folder structure improvements.
* Comment additions to code.
* Better README.md.
* `gg.py` which is a GUI for gpaw-tools is added to project. It can do all `gpawsolve.py`'s features in a graphical way!
* `gpawsolve.py` can be run solely as a command now (This is needed for a GUI project).
* All three scripts`PW-Electronic.py`, `LCAO-Electronic.py` and `PW-Optical-SingleCoreOnly.py` scripts becomes a single for-all-case script: `gpawsolve.py`.
* `PW-Electronic-changename.py` script becomes `PW-Electronic.py`.
* Spin-polarized results in `PW-Electronic-changename.py` script.
* All-electron density calculations in `PW-Electronic-changename.py`.
* CIF Export in `PW-Electronic-changename.py` script.
* Better parallel computation.
* Several XCs available for PW.
* `LCAO-Electronic.py` script.
* Strain minimization in PW only. 
* BFGS to LBFGS, Small many changes have been done.

### Preversion
* `PW-Optical-SingleCoreOnly.py` script for optical calculations.
* `PW-Electronic-changename.py` script for electronic calculations. 
* First scripts for personal usage.