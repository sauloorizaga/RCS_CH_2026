# RCS_CH_2026
A Second-Order Richardson--Convex Splitting Method for the 
Cahn-Hilliard Equation: Stability Analysis and GPU-Accelerated 
3D Computations

## Overview
We propose the Richardson--Convex Splitting (RCS) framework, 
a second-order temporal discretization based on composition 
of first-order CS operators. The scheme eliminates the parasitic 
spectral root present in BDF2-CS, requires a relaxed splitting 
parameter ($a \geq 1$ vs $a \geq 4$), and scales to $256^3$ 
3D simulations on consumer GPU hardware. At 512^3, the solver requires at least a 24 GB GPU in double precision.

## Requirements
- MATLAB with Parallel Computing Toolbox
- NVIDIA GPU (any CUDA-capable GPU)

## Contents
- `CH3D_RCSX_Solver.m` — 3D RCS solver function
- `RCS_Simulation_Script2026.m` — 3D CH simulation script (uses the above solver)
- `RCS_Comparison2D.m` — 2D comparison: RCS vs BDF2-CS
- `RCS_2D_code.m` — self contained script to run the RCS method in 2D

## Output
Running **RCS_Simulation_Script2026.m** will generate 
the 3D CH simulation and energy plot:

<img src="CH3D.png" width=400px height=400px> <img src="Energy.png" width=400px height=400px>

## Output
Running **RCS_Comparison2D.m** will generate the 2D side-by-side 
morphology comparison and energy profiles:

<img src="RCSvsBDF.png" width=500px height=400px>

<img src="RCS_BDF_Energy.png" width=500px height=400px>

## Usage
```matlab
% Parameters: N=64, Tf=10, dt=0.01, eps=0.1
% 3D simulation: N=64, Tf=10, dt=0.01, eps=0.1
% Place CH3D_RCSX_Solver.m in the same folder with RCS_Simulation_Script2026.m 
then run RCS_Simulation_Script2026.m 

% To simulate the 2D RCS vs BDF2 comparison: N=64, Tf=100, dt=0.01, eps=0.05
please run RCS_Comparison2D.m

% To run the RCS alone for 2D computations,
use the code RCS_2D_code.m  
```

## Citation
If you use this code, please cite:
Orizaga, S. (2026).
"A Second-Order Richardson--Convex Splitting Method for the 
Cahn--Hilliard Equation: Stability Analysis and GPU-Accelerated 
3D Computations"
Submitted for publication.
Code available at:
https://github.com/sauloorizaga/RCS_CH_2026

## Contact
We welcome questions, feedback, and potential collaboration 
opportunities — feel free to reach out! <br>
**Saulo Orizaga** — saulo.orizaga@nmt.edu <br>
Associate Professor of Mathematics <br>
New Mexico Institute of Mining and Technology <br>
Socorro, NM 87801, USA.
