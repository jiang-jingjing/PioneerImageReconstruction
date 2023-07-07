# Pioneer image reconstruction 
Author: Jingjing Jiang jing.jing.jiang@outlook.com

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
  
### Related publications
-  <span style="font-family:Lato; font-size:16px;">**Jiang, J.**, Di Costanzo Mata, A., Lindner, S., Charbon, E., Wolf, M. & Kalyanov, A. Dynamic time domain near-infrared optical tomog- raphy based on a SPAD camera. Biomed. Opt. Express 11, 5470 (2020). [link](https://opg.optica.org/boe/fulltext.cfm?uri=boe-11-10-5470&id=437959) </span>
- <span style="font-family:Lato; font-size:16px;">**Jiang, J.**, Di Costanzo Mata, A., Lindner, S., Zhang, C., Charbon, E., Wolf, M. & Kalyanov, A. Image reconstruction for novel time do- main near infrared optical tomography: towards clinical applications. Biomed. Opt. Express 11, 4723 (2020). [link](https://opg.optica.org/boe/fulltext.cfm?uri=boe-11-8-4723&id=433907)</span> 
[More ...](https://scholar.google.com/citations?user=hoy7VbIAAAAJ&hl=en&oi=sra) </span>
 
