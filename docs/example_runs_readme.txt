This folder contains several example runs of the model, including
single-element validation runs, cyclic strain of quarter cylinder
with symmetry boundaries to model the DMA experiments, and the example
demonstrated in the manuscript of a simple bi-layer pad hit with a
standard impactor at 3 m/s.

Folders include:
/bilayer_pad_example_fig7
	-> .inp file and run script for the impactor example shown in
	   Fig 7 of the associate manuscript.
	
/quart_cyl_cyclic
	-> Various .inp files with element definitions for
	   differing complexity meshes (single and multi-element tests),
	   designed to simulate the DMA experiment. Use for model fitting
	   and validation simulations. Also includes several python scripts
	   used to extract data from the resulting output files for use in
	   post-processing.
	   
/step_compress
	-> Includes early model development verification - a step compression
	   and hold, with varying model parameters and a summary of the
	   results (as .xlsx)
	   
/uniaxial_tension
	-> Example of uniaxial tension, used mostly for model verification
	   and development.


To run, first make sure CCX is compiles properly using the SLS+radiation
UMAT provided in this package. Make sure the paths within the run scripts
point to whereever the compiled executible resides. Next, make "run_ccx.sh"
executable (e.g., with "chmod +x run_ccx.sh" in a bash terminal), and then
simply execute "./run_ccx.sh". The results will be written to the containing
folder, and progress can be monitored using the .sta file.

Note that the codes have been tested on OpenMP parallization. 

Analyzes of the result have been conducted using either the included python
scripts in /src (for, e.g., stress traces) or by using the software package
"ccx2paraview" (from https://github.com/calculix/ccx2paraview) for
visualization within ParaView (https://www.paraview.org/). Note the results
files from these examples have not been included.