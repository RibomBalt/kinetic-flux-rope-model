# Kinetic Flux Rope Model

The model code for our paper under preparation.

## Description

Flux ropes are common magnetic structures in space plasmas with helical magnetic field lines and strong core field. Such structures can be constructed as a self-consistent solution of equilibrium Vlasov-Maxwell equations in cylindrical coordinates. For 1D scenarios, i.e. electromagnetic field and particle distribution functions only depend on cylindrical distance r, there are mature procedure to derive such solutions (e.g. [Vinogradov et al. (2016)](http://aip.scitation.org/doi/10.1063/1.4958319)). One needs to assume the particle distribution as some function of invariants of motion, and solve the Maxwell equations, i.e. the Ampere's and quasi-neutral equation with an ODE solver. 

Our model is a generalization of [Vinogradov et al., (2016)](http://aip.scitation.org/doi/10.1063/1.4958319). The original model uses electron energy, angular canonical momentum and axial momentum to construct the distribution function. We add the adiabatic invariant μ to cover the anisotropic nature of flux ropes. Note that to compute μ, the magnetic field is pre-required, so we adopted a iteration procedure to find the solution (similar to [Li et al., (2020)](http://www.nature.com/articles/s41467-020-19442-0)).

Our model takes a set of physical parameters as input. The output includes the electric field, magnetic field and electron distribution function.

## Requirements

The code is tested on MATLAB R2018a (Windows 10), but should work fine for most versions of MATLAB.

The code costs about half an hour on our computer (Intel i5-10400 CPU 2.90GHz, 16GB RAM). Parallel Computing Toolbox for MATLAB is recommended.

## Usage

The `1. params` block in `vinogradov_mur_noiteration.m` contains the parameters of this model. The meaning of these parameters are noted in Table 1 of our paper. The units are commented in the code. Many of our parameters have the same meaning as [Vinogradov et al., (2016)](http://aip.scitation.org/doi/10.1063/1.4958319).

The `inputs` variable in `2. solver & moments` block in `vinogradov_mur_noiteration.m` is the location (normalized) where the equations are solved. `inputs` are used as key points of interpolation in the subsequent iteration procedure, so it should be dense enough and cover the full area of interest. 

The `niteration` in `vinogradov_murz_int3_doiteration.m` controls the iteration times of this code. The iteration should not stop until the electromagnetic field has seldom change before and after the iteration. The default value is 3 times.

After these parameters are set, first run `vinogradov_mur_noiteration`, then `vinogradov_murz_int3_doiteration`.

The variable `final_line` is a struct that contains everything the model returns. The meaning of most useful fields are listed below:

| Field names | Meaning                                                      |
| ----------- | ------------------------------------------------------------ |
| x           | Same as "inputs", cylindrical distance normalized by `params.Lst` |
| Bz          | z component of magnetic field in nT                          |
| Bphi        | φ component of magnetic field in nT                          |
| E           | r component of electric field in mV/m                        |
| Ni / Ne     | Ion/electron number density in per cm-3                      |
| vphi        | φ component of electron velocity in km/s                     |
| vz          | z component of electron velocity in km/s                     |
| Prr/Ppp/Pzz | Pressure tensor in corresponding directions in nPa           |
| Trr/Tpp/Tzz | Tempreture tensor in corresponding directions in eV          |
| fe          | Electron distribution function. The input includes three components of velocity (in r, φ, z directions respectively, normalized by `params.v0e`) and location. |

The expected result are in `data/model_res.mat`

## License

This code is covered under Apache License 2.0.
