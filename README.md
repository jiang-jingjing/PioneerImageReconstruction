# Pioneer image reconstruction 
## Introduction

The  is a package based on MATLAB for the time-domain 
near infrared optical tomography system Pioneer.

It incorporates various image reconstruction methods, forward models and
 tissue geometries, source detector arrangements, etc.

### Preprocessing
- Step 1: run [timing_data_correction_par.m](https://github.com/jiang-jingjing/PioneerImageReconstruction/blob/master/Piccolo_data_correction/timing_data_correction_z3_par.m)
to get corrected timing responses in .mat format
- Step 2: further preprocessing for the target forward models
### Forward models
Methods for simulating light propogation or diffusion inside turbid media.
#### Numerical models include:
 
  - Monte Carlo Methods
  - Finite Element Methods
  
  
