# MATLAB MPC Tutorials using CasADi

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

## MPC formulations

This code assumes a linear, state-space model generated from the plasma jet testbed. The two inputs $u$ (manipulated variables) are power ($P$, in Watts) and carrier gas flow rate ($q$, in standard liters per minute). The two measured outputs $y$ (control variables, assuming state feedback $y=x$) are maximum surface temperature ($T$, in degrees Celsius) and total optical emission intensity ($I$, in arbitrary units).

$$x = [T, I]^\top$$

$$u = [P, q]^\top$$

$$x^+ = Ax + Bu$$

where $x^+$ denotes the successor state and $A$ and $B$ are the model matrices obtained through system identification via the `n4sid` function in MATLAB.

### Nominal (Tracking) MPC
The objective of this controlled system is to track reference values of the temperature, i.e., over a period of time, we wish to control our system to attain desired temperature set points.

The mathematical formulation of the optimal control problem is

```math
\begin{align}
\min_{\mathbf{x},\mathbf{u}} &~~ \|y_{N}-y_{ref}\|_{Q_N} + \sum_{k=0}^{N-1} \|y_{k}-y_{ref}\|_{Q}, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ y_k = x_k^1, \\
&~~ x_0 = x(0), \\
&~~ [25^\circ\text{C}, 0 \text{ arb. units}]^\top \le x_k \le [45^\circ\text{C}, 80 \text{ arb. units}]^\top, \\
&~~ [1.5 \text{ W}, 1.5 \text{ SLM}]^\top \le x_k \le [5.0 \text{ W}, 5.0 \text{ SLM}]^\top, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```

**Note:** The current implementation only tracks one of the states, so $y_k = x_k^1 = T_k$.

### Offset-free MPC
**TODO: edit** 
The objective of this controlled system is to track reference values of the temperature, i.e., over a period of time, we wish to control our system to attain desired temperature set points.

The mathematical formulation of the optimal control problem is

```math
\begin{align}
\min_{\mathbf{x},\mathbf{u}} &~~ \|y_{N}-y_{ref}\|_{Q_N} + \sum_{k=0}^{N-1} \|y_{k}-y_{ref}\|_{Q}, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ y_k = x_k^1, \\
&~~ x_0 = x(0), \\
&~~ [25^\circ\text{C}, 0 \text{ arb. units}]^\top \le x_k \le [45^\circ\text{C}, 80 \text{ arb. units}]^\top, \\
&~~ [1.5 \text{ W}, 1.5 \text{ SLM}]^\top \le x_k \le [5.0 \text{ W}, 5.0 \text{ SLM}]^\top, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```

**Note:** The current implementation only tracks one of the states, so $y_k = x_k^1 = T_k$.

### Economic MPC
**TODO: edit** 
The objective of this controlled system is to track reference values of the temperature, i.e., over a period of time, we wish to control our system to attain desired temperature set points.

The mathematical formulation of the optimal control problem is

```math
\begin{align}
\min_{\mathbf{x},\mathbf{u}} &~~ \|y_{N}-y_{ref}\|_{Q_N} + \sum_{k=0}^{N-1} \|y_{k}-y_{ref}\|_{Q}, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ y_k = x_k^1, \\
&~~ x_0 = x(0), \\
&~~ [25^\circ\text{C}, 0 \text{ arb. units}]^\top \le x_k \le [45^\circ\text{C}, 80 \text{ arb. units}]^\top, \\
&~~ [1.5 \text{ W}, 1.5 \text{ SLM}]^\top \le x_k \le [5.0 \text{ W}, 5.0 \text{ SLM}]^\top, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```

### Multistage MPC
**TODO: edit** 
The objective of this controlled system is to track reference values of the temperature, i.e., over a period of time, we wish to control our system to attain desired temperature set points.

The mathematical formulation of the optimal control problem is

```math
\begin{align}
\min_{\mathbf{x},\mathbf{u}} &~~ \|y_{N}-y_{ref}\|_{Q_N} + \sum_{k=0}^{N-1} \|y_{k}-y_{ref}\|_{Q}, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ y_k = x_k^1, \\
&~~ x_0 = x(0), \\
&~~ [25^\circ\text{C}, 0 \text{ arb. units}]^\top \le x_k \le [45^\circ\text{C}, 80 \text{ arb. units}]^\top, \\
&~~ [1.5 \text{ W}, 1.5 \text{ SLM}]^\top \le x_k \le [5.0 \text{ W}, 5.0 \text{ SLM}]^\top, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```


### Minimum Time MPC
**TODO: edit** 
The objective of this controlled system is to track reference values of the temperature, i.e., over a period of time, we wish to control our system to attain desired temperature set points.

The mathematical formulation of the optimal control problem is

```math
\begin{align}
\min_{\mathbf{x},\mathbf{u}} &~~ \|y_{N}-y_{ref}\|_{Q_N} + \sum_{k=0}^{N-1} \|y_{k}-y_{ref}\|_{Q}, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ y_k = x_k^1, \\
&~~ x_0 = x(0), \\
&~~ [25^\circ\text{C}, 0 \text{ arb. units}]^\top \le x_k \le [45^\circ\text{C}, 80 \text{ arb. units}]^\top, \\
&~~ [1.5 \text{ W}, 1.5 \text{ SLM}]^\top \le x_k \le [5.0 \text{ W}, 5.0 \text{ SLM}]^\top, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```
