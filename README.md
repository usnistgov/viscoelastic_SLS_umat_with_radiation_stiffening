# Viscoelastic Standard Linear Solid User Material with Thermal and Radiation Dose Effects

## Purpose
This repository contains code and calibration data for a linear viscoelastic (Standard Linear Solid rheological model) finite element user material that incorporates the effects of temperature and radiation dosage via empirical fitting to measured data.

The model is a research-grade code and while it should be relatively plug-and-play for those comfortable with deploying finite element user materials in Calculix/CrunchiX or Abaqus, is it not standalone software and requires operator skill to use successfully. 

## Content
Since this is not a standalone package the content and workflow is somewhat different from a more conventional software release, and the content of each subdirectory is described below.

- The top level directory of this project contains several documentary files, including this readme, license information, and other metadata.
- The `bin` directory contains an example complied binary executable for CalculiX/CrunchiX that include the irradiated foam  user material. Warning: this may not execute on in your environment (built under Ubuntu 20.04.6 LTS on Windows Subsystem for Linux 2 x64).
- The `calibration_data` directory includes the curated DMA data used as input to fit the model, as well as the calibration data handling, fitting, metadata, and scripts.
- The `docs` directory contains text-based documentation for compiling from source code and running the example cases. It also includes note on the theory used to implement the viscoelasticity model in a `.pdf`, and the spreadsheet (`.xlsx`) file with validation data, as shown in Fig. A8 of the manuscript.
- The `example_runs` directory has subfolder for four example test cases using the user material under various conditions. Specifically, the example padding application shown in Fig 7 of the manuscript, cyclic compression, stepwise compression, and uniaxial tension.
- The `src` directory includes source code (.f and .py) for the material model, data handling, and results analysis. 

The Calculix source code (`umat.f`) and input file (`*.inp`) are syntactically similar to the source code used for Abaqus - it is highly plausible to port this code to Abaqus, likely with relatively few changes, although we have not test it (yet).   

## Running the user material and example models
See the `docs` and `example_runs` subdirectories for more complete directions. You may also wish to consult CalculiX documentation (see <http://www.dhondt.de/>).

Briefly, to run an example you will need to:
1. Compile CalculiX from source adding the user material to the compiled executable. Note: you may be able to skip this step by using the pre-compiled binary in `bin` if your runtime environment is suitably similar to ours.
2. Update the run script (typically, `./run_ccx.sh`) for the example with the correct path to the CCX executable.
3. Run the `.sh` file in a bash terminal.
4. Results can be processed and analyze with included Python scripts, converted to be visualized ParaView ("ccx2paraview", <https://www.paraview.org/>), or processed with your own tools.

## Other links
For the data used in the development and calibration of the model see: <https://doi.org/10.18434/mds2-2989>.

For more details and description of the experiment see: <https://doi:10.1016/j.matdes.2023.112381>.

Please cite as:
> Landauer, A. K. et al. Unintended consequences: Assessing thermo-mechanical changes in vinyl nitrile foam due to micro-computed X-ray tomographic imaging. Materials & Design 112381 (2023) doi:10.1016/j.matdes.2023.112381.


## Contact and support
For questions, please open a new entry in the "Issues" tab. If needed, you can also find authors' contact information via the associated paper (see above). 

The corresponding author is Alex Landauer (NIST MML Materials Measurement Science Division, Security Technologies Group). Orion Kafka and Newell Moser (NIST MML Applied Chemicals and Materials Division) are the primary developers of the code.

---

