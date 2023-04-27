# MPC_tutorials

Tutorials on various model predictive control (MPC) implementations in MATLAB. Created for the Mesbah Lab at the University of California, Berkeley.

The main scripts to run the code are the files `run_xxx.m` where the `xxx` indicates a different formulation of the MPC problem.

This repository includes MPC formulations of:
* nominal MPC - no disturbances with simple quadratic reference tracking cost
* offset-free MPC - quadratic reference tracking with offset correction
* economic MPC - "economic" cost, in which the objective is to achieve a target cumulative thermal dose delivery; the cost is formulated w.r.t. tracking a desired dose
* multistage MPC - a robust formulation of the economic MPC, in which uncertainty is incorporated via scenario tree and the worst-case bounds are used to generate the scenarios
* minimum time MPC - minimum time cost in which the thermal dose delivery is formulated such that the treatment time is directly minimized vs the tracking formulation in the economic and multistage formulations


MATLAB Livescripts may be found in the `Live Scripts` folder to provide an "interactive" demo of the code. If you use these scripts, you should open them and then make sure your working directory is the main directory where the `run_xxx.m` files are.

These scripts are created with a project-oriented view. Thus, the controller, closed-loop simulation, etc. codes are abstracted away to files located under the `utils` folder. To create new controller formulations, you should modify/create your own version of the controller scripts under the `utils` folder. Furthermore, the `config` folder contains scripts that contain problem-specific configurations. As such, those can be used as a template/starting-point when transferring your codes from one system to another.

Example output from the scripts are provided in `example_results.pptx`.
