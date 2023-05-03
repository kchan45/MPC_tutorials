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

### Notation
A controlled system may run indefinitely, but for simulation purposes, we create a finite simulation horizon $N_{sim}$.

For evolution of the states within the simulation horizon (aka as a part of the real trajectory of the true system), the state at time $t$ is denoted as $x(t)$. The inputs and outputs are defined similarly: the input at time $t$ is $u(t)$; the output at time $t$ is $y(t)$.

The state evolution dynamic representation of the true system can be written as
```math
x(t+\Delta t) = f_p(x(t), u(t)),
```
and the output representation of the true system can be written as
```math
y(t) = h_p(x(t), u(t)),
```
where no assumptions are made on noise, disturbances, uncertainty, etc.

MPC solves an optimal control problem (OCP) at every time instant of the simulation horizon in order to obtain an optimal trajectory of inputs. Thus, it utilizes a model to predict the future states and outputs (as a result of applying certain inputs). Often, the prediction is made up until some user-defined prediction horizon $N\le N_{sim}$.

For the prediction of the states within an OCP, the state evolution model is denoted as
```math
x_{t+\Delta t} = f(x_{t}, u_{t}),
```
and the output model is denoted with
```math
y_{t} = h(x_{t}, u_{t}),
```
where a prediction is denoted with a subscripted time instant.

The following are equivalent exemplify the connection between the two notations
```math
x_{t_1} = x(t_1|t),
```
each of which designates the predicted state at time $t_1$ in the prediction horizon $N$ given the computation of the OCP at time $t$ in the simulation horizon $N_{sim}$.

In multi-state/input/output systems, the $n$-th state is denoted by a superscript, i.e., the first, third, and fifth state at time $t$ in the simulation horizon are denoted by $x^1(t)$, $x^3(t)$, and $x^5(t)$, respectively, while the first, third, and fifth state *prediction* at time $t$ of the prediction horizon are denoted by $x_t^1$, $x_t^3$, and $x_t^5$, respectively.

### System/Problem Information
This code assumes a linear, state-space model generated from the plasma jet testbed. The two inputs $u$ (manipulated variables) are power ($P$, in Watts) and carrier gas flow rate ($q$, in standard liters per minute). The two measured outputs $y$ (control variables, assuming state feedback $y=x$) are maximum surface temperature ($T$, in degrees Celsius) and total optical emission intensity ($I$, in arbitrary units).

$$x = [T, I]^\top$$

$$u = [P, q]^\top$$

$$x^+ = Ax + Bu$$

where $x^+$ denotes the successor state and $A$ and $B$ are the model matrices obtained through system identification via the `n4sid` function in MATLAB.

The states are constrained as

$$x_{\min} = [25^\circ\text{C}, 0 \text{ arb. units}]^\top$$

$$x_{\max} = [45^\circ\text{C}, 80 \text{ arb. units}]^\top$$

The inputs are constrained as

$$u_{\min} = [1.5 \text{ W}, 1.5 \text{ SLM}]^\top$$

$$u_{\max} = [5.0 \text{ W}, 5.0 \text{ SLM}]^\top$$


### Nominal (Tracking) MPC
The objective of this controlled system is to track reference values of the temperature, i.e., over a period of time, we wish to control our system to attain desired temperature set points.

The mathematical formulation of the optimal control problem is

```math
\begin{align}
\min_{\mathbf{x},\mathbf{u}} &~~ \|y_{N}-y_{ref}\|_{Q_N}^2 + \sum_{k=0}^{N-1} \|y_{k}-y_{ref}\|_{Q}^2, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ y_k = x_k^1, \\
&~~ x_0 = x(0), \\
&~~ x_{\min}^\top \le x_k \le x_{\max}, \\
&~~ u_{\min} \le u_k \le u_{\max}, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```

**Note:** The current implementation only tracks one of the states, so $y_k = x_k^1 = T_k$.

### Offset-free MPC
The objective of this controlled system is to track reference values of the temperature *without offset*, i.e., over a period of time, we wish to control our system to attain desired temperature set points *without a gap between our achieved steady state versus the reference*.

There are two well-established strategies to address offset-free tracking:
1. Rawlings(and Pannocchia)-style formulation [1](https://doi.org/10.1002/aic.690490213), [2](https://doi.org/10.1109/ECC.2015.7330597), [3](https://doi.org/10.1016/j.ifacol.2015.11.304)
2. Limon-style formulation [4](https://doi.org/10.1007/978-3-642-01094-1_26), [5](https://doi.org/10.1016/j.automatica.2008.01.023), [6](https://doi.org/10.1109/CDC.2009.5400618)

The tutorial code currently implements the Limon-style, and the Rawlings-style may be added in a future update of this code. The mathematical formulation of the optimal control problem for the Limon-style is

```math
\begin{align}
\min_{\mathbf{u},\theta} &~~ \|x_{N}-x_{s}\|_{Q_N}^2 + \|y_{N}-y_{s}\|_{P_N}^2 + \sum_{k=0}^{N-1} \|x_{k}-x_{s}\|_{Q} + \|u_{k}-u_{s}\|_{R}, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ x_0 = x(0), \\
&~~ \theta = (x_{s}, u_{s}, y_{s}), \\
&~~ x_{s} = g_{x}(\theta) = f(x_{s},u_{s}) = Ax_{s} + Bu_{s}, \\
&~~ y_{s} = g_{y}(\theta) = h(x_{s},u_{s}) = Cx_{s} + Du_{s}, \\
&~~ [25^\circ\text{C}, 0 \text{ arb. units}]^\top \le x_k \le [45^\circ\text{C}, 80 \text{ arb. units}]^\top, \\
&~~ [1.5 \text{ W}, 1.5 \text{ SLM}]^\top \le x_k \le [5.0 \text{ W}, 5.0 \text{ SLM}]^\top, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```

**Note:** The current implementation only tracks one of the states, so $y_k = x_k^1 = T_k$.

### Economic MPC
The objective of this controlled system is to attain a desired "thermal dose" of plasma. The "thermal dose" quantifies the amount of thermal effects delivered to a surface. We use the *cumulative equivalent minutes (CEM)* as the metric to quantify this thermal dose. The CEM is derived from hypothermia treatments and is defined as

```math
\mathrm{CEM}(t) = \int_{0}^{t} K^{(T_{\mathrm{ref}} - T(\tau))} d\tau,
```
where $K$ is a temperature-dependent exponential base that is related to heat transfer properties of a substrate/surface with the plasma. It is relevant to understand that the CEM is *cumulative* (non-decreasing).

The mathematical formulation of the optimal control problem is

```math
\begin{align}
\min_{\mathbf{x},\mathbf{u}} &~~ \|\mathrm{CEM}_{N}-\mathrm{CEM}_{sp}\|^2, \\
\text{s.t.} &~~ x_{k+1} = Ax_k + Bu_k, \\
&~~ x_0 = x(0), \\
&~~ [25^\circ\text{C}, 0 \text{ arb. units}]^\top \le x_k \le [45^\circ\text{C}, 80 \text{ arb. units}]^\top, \\
&~~ [1.5 \text{ W}, 1.5 \text{ SLM}]^\top \le x_k \le [5.0 \text{ W}, 5.0 \text{ SLM}]^\top, \\
&~~ \mathbf{x} = \{x_k\}_{k=1}^{N}, ~~\mathbf{u} = \{u_k\}_{k=0}^{N-1}, \\
&\quad\quad \forall k = 0, \ldots, N-1,
\end{align}
```

### Multistage MPC
The objective of this controlled system is also to attain a desired "thermal dose" of plasma (see Economic MPC section for details on the definition of the thermal dose). However, multistage (or scenario-based) MPC is considered a "robust" MPC, which explicitly considers uncertainty quantification in the formulation. Multistage MPC considers propagates uncertainty via a scenario tree []().

*Note:* This is not the only robust MPC formulation. An alternative robust MPC formulation (and more common/initial version of robust MPC) is **tube MPC**, which propagates uncertainty of the system as a tube of trajectories.

**TODO: edit**
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
