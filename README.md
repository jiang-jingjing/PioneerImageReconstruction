# Pioneer image reconstruction 
created by Jingjing Jiang jing.jing.jiang@outlook.com

## Introduction

The  is a package based on MATLAB for the time-domain 
near infrared optical tomography system Pioneer.

It incorporates various image reconstruction methods, forward models and
 tissue geometries, source detector arrangements, etc.

The code was implemented in Linux OS. If you are using Windows OS, you may need to adjust a few things, e.g., \ or / used for defining paths.
### Preprocessing
- Step 1: run [timing_data_correction_par.m](https://github.com/jiang-jingjing/PioneerImageReconstruction/blob/master/Piccolo_data_correction/timing_data_correction_z3_par.m) 
to get corrected timing responses in .mat format. (Requirement: Parallel computing toolbox installed in MATLAB, otherwise run [timing_data_correction_z3.m](https://github.com/jiang-jingjing/PioneerImageReconstruction/blob/master/Piccolo_data_correction/timing_data_correction_z3.m) ).  
- Step 2: further preprocessing for the target forward models
### Forward models
Methods for simulating light propogation or diffusion inside turbid media.
#### Numerical models include:
 
  - Monte Carlo Methods
  see the [example](https://github.com/jiang-jingjing/PioneerImageReconstruction/blob/master/exampleMC.m)
  - Finite Element Methods
  see the [example](https://github.com/jiang-jingjing/PioneerImageReconstruction/blob/master/exampleFEM.m)
  
