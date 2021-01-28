# modelSSNMR
Solid-State NMR restraints for molecular dynamics (GROMACS patch)

- [General remarks](#general-remarks)
- [ssNMR models and theory](#ssnmr-models-and-theory)
- [Installation](#installation)
- [Running ssNMR restrained simulations](#running-ssnmr-restrained-simulations)

## General remarks

The modelSSNMR patch modifies a series of C source and include files in the [GROMACS](https://www.gromacs.org/) suite (version 4.6.7) in order to run restrained Molecular Dynamics (MD) simulations using data from solid-state nuclear magnetic resonance (ssNMR).

If using or referencing this methodology, please cite:

*Sanz-Hernández, M. et al. Accurate Determination of Conformational Transitions in Oligomeric Membrane Proteins. Sci. Rep. 6, 23063;* doi: [10.1038/srep23063](https://www.doi.org/10.1038/srep23063) (2016)

*De Simone, A. et al. Structural Dynamics and Conformational Equilibria of SERCA Regulatory Proteins in Membranes by Solid-State NMR Restrained Simulations. Biophys. J. 106(12): 2566–76;* doi: [10.1016/j.bpj.2014.03.026](https://www.doi.org/10.1016/j.bpj.2014.03.026) (2014)

## ssNMR models and theory

The ssNMR data used here to restrain MD simulations are [chemical shift anisotropy](http://triton.iqfr.csic.es/guide/eNMR/proteins/CSA.html) (CSA) and [dipolar couplings](https://en.wikipedia.org/wiki/Magnetic_dipole%E2%80%93dipole_interaction) (DC). These data provide very powerful topological information regarding the orientation of chemical groups in proteins, and are particularly useful in the context of membrane proteins.

### Definition of CSA and DC models

In the protein backbone, the CSA and DC can be calculated given the 3D coordinates of atoms, according to the following structural models:

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/MD_SSNMR_Calculator/SSNMR_models.png" width="800" title="SSNMR_Models">

The CSA is modelled as a rank-2 tensor centered on the amide nitrogen atom, with the following equation:

![equation0](https://latex.codecogs.com/svg.latex?\delta_{15_N}=\delta_{11}\times\sin^2(\alpha-17)\times\sin^2\beta&plus;\delta_{22}\times\cos^2\beta&plus;\delta_{22}\times\cos^2(\alpha-17)\times\sin^2\beta)

where *&delta;<sub>11</sub>*, *&delta;<sub>22</sub>* and *&delta;<sub>33</sub>* are respectively set to 64.0, 76.0, 216.9 ppm for non-glycine residues and 46.5, 66.3, 211.6 ppm for glycine. *&alpha;* and *&beta;* are the Euler angles (in degrees) used to transform from the laboratory frame to the principal axis frame.

The DC is only dependent on the length of the covalent bond and its angle &theta; with respect to the external magnetic
field B<sub>0</sub>:

![equation1](https://latex.codecogs.com/svg.latex?DC_{NH}=\frac{1}{2}\zeta_{DC}(3\cos^2\theta-1))

where *&zeta;<sub>DC</sub>* is set to 10.52 kHz for the <sup>15</sup>N-<sup>1</sup>H amide bond and -22.68 kHz for the <sup>13</sup>C<sub>&alpha;</sub>-<sup>1</sup>H<sub>&alpha;</sub> bond.

### Restraining methodology

The objective of the restrained sampling algorithm is to maximize the agreement between experimental and back-calculated CSA, DC data from the simulation.

In order to do so, we introduce a harmonic potential (V<sub>ssNMR</sub>) that is added to the molecular dynamics force field:

![equation2](https://latex.codecogs.com/svg.latex?V_{ssNMR}=k\sum_i(Obs^{Calc}_i-Obs^{Exp}_i)^2)

where *i* runs over all the available data instances. *Obs* refers to any observable (CSA and DC in our case), with *Obs<sup>Exp</sup>* referring to its experimental value and *Obs<sup>Calc</sup>* to the back-calculated. *k* is a system-dependent harmonic force constant that determines the strength of the restraint. The potential is normally introduced with a flat bottom, whereby no force is applied if the difference between *Obs<sup>Exp</sup>* and *Obs<sup>Calc</sup>* is below an error value *&epsilon;*.

In addition, the algorithm must consider that the observables are the result of [ensemble averaging](https://en.wikipedia.org/wiki/Ensemble_average_(statistical_mechanics)) across all the conformations that compose the experimental sample.

To account for this, we perform **replica-averaged restrained MD**, where multiple copies (known as replicas) of the system are simulated simultaneously and independently. At every step, the observables are back-calculated for each replica and averaged amongst the *M* replicas:

![equation3](https://latex.codecogs.com/svg.latex?Obs^{Calc}_i=\frac{1}{M}\sum_{m}^{M}Obs^{Calc}_{i,m})

There is an additional level of averaging possible. In homo-oligomeric proteins, there are  multiple copies of the same molecule interacting, which further contribute to the ensemble averaging of the experimental data. To account for this property, we implemented internal averaging amongst the P protomers that make up the quaternary protein structure, and combined this with the replica averaging scheme:

![equation4](https://latex.codecogs.com/svg.latex?Obs^{Calc}_i=\frac{1}{M}\frac{1}{P}\sum_{m}^{M}\sum_{p}^{P}Obs^{Calc}_{i,m,p})

### Restraining forces

At every step of the simulation, after calculating the value of the restraining potential, *V<sub>ssNMR</sub>*, a resulting force (***F<sub>ssNMR</sub>***) is applied to maximize the agreement with experimental data:

![equation5](https://latex.codecogs.com/svg.latex?\vec{F}_{ssNMR}=-\frac{dV_{ssNMR}}{d\vec{r}})

where ***r*** are the atomic coordinates of each atom involved in the restraint.

Applying the chain rule, we can further separate the two contributions to the force:

![equation6](https://latex.codecogs.com/svg.latex?\vec{F}_{ssNMR}=-\frac{dV_{ssNMR}}{dObs^{Calc}}\frac{dObs^{Calc}}{d\vec{r}})

## Installation

The modelSSNMR patch is applied to an unmodified version of **GROMACS 4.6.7**, which can be downloaded from the [GROMACS site](http://www.gromacs.org/Downloads_of_outdated_releases).

In order to patch your GROMACS distribution, download the modelSSNMR patch and navigate to its directory. From there, you may run:

```
./patch-gromacs-4.6.7 [GROMACS_src_dir]
```

replacing ```[GROMACS_src_dir]``` with the directory where your uncompressed GROMACS source files are.

This will apply all necessary changes to GROMACS, which can then be installed following the [GROMACS installation instructions](http://www.gromacs.org/Documentation_of_outdated_versions/Installation_Instructions_4.6).'

In order to run **replica-averaged restrained MD** (with the ```-multi``` flag), GROMACS-modelSSNMR must be installed with MPI enabled.

As a result of the patch, the original ```[distance_restraints]``` module of GROMACS is turned off. However, we provide another option for distance restraints (see below).

## Running ssNMR restrained simulations

The preparation of the system is identical to a general protocol for GROMACS MD simulations.

Experimental data are instructed via the file **ssnmr**, which should be located in the simulation working directory

### ssnmr file format

Format of the **ssnmr** file
Each line should include one experimental data with the following elements separated by whitespace.

For **CSA**:
- Identifier
- First_atom
- Second_atom
- Third_atom
- Experimental_value (ppm)
- Error_tolerance (ppm, optional – default: 5.0)
- Tensor_magnitude_11 (ppm, optional – default: 64.0)
- Tensor_magnitude_22 (ppm, optional – default: 76.0)
- Tensor_magnitude_33 (ppm, optional – default: 216.9)

For **DC**:
- Identifier
- First_atom
- Second_atom
- Experimental_value (kHz)
- Error_tolerance (kHz, optional – default: 0.5)
- Dipolar coupling constant (optional – default: 10.52 kHz, for amide N-H)

For **Distances**:
- Identifier
- First_atom
- Second_atom
- Experimental_value (nm)
- Error_tolerance (nm, optional – default: 0.05)

Identifiers are as follows: 
- 0 = CSA
- 1 = DC
- 2 = Distance

More specifically, each line in the ssnmr file should represent 1 experimental datum. The lines are different if the data is CSA or DC:

In the case of CSA (identifier 0), for a given peptide plane of residue *i* 
- First_atom = N (*i*)
- Second_atom = HN (*i*)
- Third_atom = CO (*i-1*)

The error tolerance can optionally be modified, as well as the magnitude of the CSA tensor (e.g. for Glycine residues)

Example of experimental data for a set of <sup>15</sup>N CSA 

```
# CSA data without specified errors (default error of 5.0)
0 121 122 119 78.000
0 67 68 65 71.000

# CSA data with custom errors (optional)
0 67 68 65 71.000 7.5
0 100 101 98 79.000 7.5

# CSA data with errors and modified tensor magnitudes (optional, e.g. for Glycine)
0 140 141 138 79.000 5.0 46.5 66.3 211.6
0 154 155 152 77.000 7.5 46.5 66.3 211.6
```

In the case of DC (identifier 1), for a given peptide plane of the residue *i* 
- First_atom = N (*i*)
- Second_atom = HN (*i*)

Example of experimental data for a set of <sup>15</sup>N-<sup>1</sup>H DC:


```
# DC data without specified errors (default error of 0.5)
1 67 68 4.900
1 121 122 4.600

# DC data with specified errors (optional)
1 140 141 5.100 1.0
1 154 155 3.500 1.0

# DC data with specified errors and kDC constant (optional, e.g. for 13Cα-1Hα)
1 178 179 3.500 0.5 -22.68
1 189 190 5.100 1.0 -22.68
```

In the case of distances (identifier 2), for any two given atoms of interest:

- First_atom = Atom 1
- Second_atom = Atom 2

Example of experimental data for a set of distances (expressed in nm)

```
1 62 210 0.49
1 401 152 1.50
1 140 101 0.49
```

### mdp file options

The input mdp file is identical to the standard GROMACS, with one additional section corresponding to the ssNMR restraints:

```

; SSNMR restraints - CSA and DC - by Maximo Sanz-Hernandez & Alfonso De Simone

; For reference see Scientific Reports 6, 23063 (2016) https://doi.org/10.1038/srep23063

; For imposing the restraint please select
; SSNMR restraints type: No, Simple, Anneal
SSNMRre = Anneal
; CSA restraints Force constant
CSAre-fc = 50
; DC restraints Force constant
DCre-fc = 2000
; CA-HA DC restraints Force constant
hDCre-fc = 0
; Distance restraints Force constant
DISre-fc = 10
; Replica averaging: between different simulations in parallel (run with the -multi flag in mdrun)
; Do replica averaging (for CSA, DC, hDC , DIS) : No or Yes
CSAre-replica_aver = yes
DCre-replica_aver = yes
hDCre-replica_aver = no
DISre-replica_aver = no
; Monomer averaging: between monomer units of an oligomeric protein
; Do monomer averaging (for CSA, DC, hDC , DIS) : No or Yes
CSAre-monomer_aver = yes
DCre-monomer_aver = yes
hDCre-monomer_aver = no
DISre-monomer_aver = no
; Option to provide signed or absolute values of the DC (Absolute val. by default)
; Back-calculate the signed value of DC and hDC : No or Yes
DCre-signed = no
hDCre-signed = no
; Periodicity: the number of atoms in each monomer
monomer_periodicity = 895
; Number of monomer units per oligomer
monomer_units = 5
; SSNMR Annealing. It operates in case Anneal is selected for SSNMR ; SSNMR annealing points
SSNMR_anneal_points = 2
; SSNMR annealing times
SSNMR_anneal_times = 0 20000
; CSA annealing forces
CSA_anneal_forces = 0 50
; DC annealing forces
DC_anneal_forces = 0 2000
; hDC annealing forces
hDC_anneal_forces = 0 0
; DIS annealing forces
DIS_anneal_forces = 0 10

```

The switch of the restrain is ```SSNMRre```.

Case ```No```, the restrain is turned off.

Case ```Simple```, the restrain runs with constant forces for CSA, DC and DIS. These are ```CSAre-fc```, ```DCre-fc``` and ```DISre-fc``` respectively. 

Case ```Anneal```, the restraint runs with variable forces along the simulation. This is the most usual type as often you would like the forces to be modular in the trajectory. For example, if you are running cycles of annealing, or if you are gently imposing the force at the beginning of the simulation. The use of the ```Anneal``` option is straightforward.

In the example provided above we have 2 points. The two times are 0 and 20000 ps. These are the endpoints of the linear interpolation of the force. In the case of CSA the force will start from 0 and raise up to 50 within the first 20000 ps. After the second time, the force will be constant at 50.

Similarly in the case of DC the force (in the example) will raise form 0 (at time 0) to 2000 (at time 20000 ps).

Replica averaging is switched on/off with the ```replica_aver``` tag. It must be specified for each restraint type separately (so some restraints can be averaged and some not).

Monomer (or internal) averaging is switched on/off with the ```monomer_aver``` keyword. If switched on, the experimental data is averaged between different copies of identical molecules in an oligomeric protein (within each simulation). You must specify the number of atoms of a single protomer with ```monomer_periodicity``` and the number of protomers in total (e.g. 3 for a trimer) with ```monomer_units```. The atom numbers in the ```ssnmr``` file must correspond to the first copy of the molecule, and all the units must be identical contiguous in the structure file (.gro).

Example:

ssnmr file:

```
1 67 68 4.900
```

mdp file: 

```
DCre-monomer_aver = yes
monomer_periodicity = 895
monomer_units = 5
```

Result:

The DC will be calculated five times, between atoms:

67-68

962-963

1857-1858

2752-2753

3647-3648

and averaged.

Replica averaging is applied among every running parallel simulations, by using the ```–multi``` flag in mdrun. The experimental data will be averaged across all the independent replicas. 

The dipolar coupling can be restrained as the absolute value (by default) or including its sign (useful if the sign is known), which is activated with the ```DCre-signed``` option

All the ```hDC``` options can be ignored in the current version - this is an experimental feature.

### Running a restrained simulation

To run a restrained simulation follow the usual GROMACS pipeline:

**grompp** step:

```
grompp_ssnmr –c conf.gro –f md_restr.mdp
```

**mdrun** step (make sure there is a **ssnmr** file in the working directory):

```
mdrun_ssnmr -v
```

if you want to average across 4 parallel replicas:

```
mpirun –np 4 mdrun_ssnmr -v -multi 4
```

```mdrun_ssnmr``` will write a file called **Rho** containing the Q-factors (agreement between experiments and simulation) for each data type (except distances).

```
       0.000 Q_CSA = 0.2543; Q_hDC = 0.0000; Q_DC = 0.3323
       1.000 Q_CSA = 0.2666; Q_hDC = 0.0000; Q_DC = 0.8542
       2.000 Q_CSA = 0.2550; Q_hDC = 0.0000; Q_DC = 0.3748
       3.000 Q_CSA = 0.2612; Q_hDC = 0.0000; Q_DC = 0.6519
       4.000 Q_CSA = 0.2428; Q_hDC = 0.0000; Q_DC = 0.6153
       5.000 Q_CSA = 0.2765; Q_hDC = 0.0000; Q_DC = 0.2477
       6.000 Q_CSA = 0.2445; Q_hDC = 0.0000; Q_DC = 0.2614
       7.000 Q_CSA = 0.2666; Q_hDC = 0.0000; Q_DC = 0.4191
       8.000 Q_CSA = 0.2491; Q_hDC = 0.0000; Q_DC = 0.3154
       9.000 Q_CSA = 0.2682; Q_hDC = 0.0000; Q_DC = 0.3486
```

If you run across multiple CPUs, domain decomposition must be used (this is faster and the default anyway), restraints do not work with particle decomposition. GPU acceleration is compatible with the restraints.
