# Runge-Kutta Method Implementation

This repository contains code that implements the Runge-Kutta methods proposed in the paper  
**"Semi-implicit-explicit Runge-Kutta Method for Nonlinear Differential Equations."**  
https://arxiv.org/pdf/2504.09969  
It also includes scripts to reproduce several tables from the paper.

## Repository Structure

### `RKmethod/`

Contains implementations of various Runge-Kutta methods, including the semi-implicit-explicit schemes proposed in the paper.

| File | Purpose |
|------|---------|
| `SemiIMEXRungeKutta.m` | Main semi-IMEX solver; selects Butcher tableau by `[order, stage, index]` and integrates $y' = f_{\text{ex}} + G(t,y)\,y$ |
| `semIIMEXRungeKutta.m` | Variant entry point for the same semi-IMEX method |
| `semiIMEXRungeKuttaTrapezoid.m` | Trapezoid-rule semi-IMEX method (2-stage, 2nd-order) |
| `RKOneStep.m` | Single time-step kernel used by the semi-IMEX solver |
| `RKOneStepIMEX2.m` | Specialised single-step for IMEX-RK2, reusing the pre-factored matrix $L$ |
| `IMEXRungeKutta2.m` | 2nd-order IMEX RK for $y' = f_{\text{ex}} + g_{\text{im}}$ with time-independent $G$ |
| `IMEXRungeKutta3.m` | 3rd-order IMEX RK for the same split form |
| `RungeKutta4.m` | Classical 4th-order explicit Runge-Kutta, used as a reference solver |

### `test/`

Includes scripts to perform convergence tests for different Runge-Kutta methods, also showing the simplest example of using the function in RKmethod/.

| File | Purpose |
|------|---------|
| `test_semiIMEXRungeKutta.m` | Convergence test for the semi-IMEX methods on the scalar ODE |

### `tools/`

Contains utility functions used in the tests, such as:

- Construction of differential matrices using the finite difference method
- Conversion of numerical results into LaTeX-formatted tables

| File | Purpose |
|------|---------|
| `diffmat.m` | Builds finite-difference differentiation matrices $[D_0, \ldots, D_m]$ of any order |
| `to_latex_convergence_table.m` | Formats an error array into a LaTeX `tabular` convergence table |

### `paper/`

Includes scripts to reproduce Tables 11, 12, and 14 from the paper.

| File | Purpose |
|------|---------|
| `paper_scalar_equation.m` | Reproduces Table 11: convergence study for the scalar ODE $y' = \cos(t)\,y + (-y + \cos(t))\,y$ |
| `paper_nonlinear_diffusion.m` | Reproduces Table 12: convergence study for the nonlinear diffusion PDE with periodic BCs |
| `paper_Cahn_Hilliard.m` | Reproduces Table 14: convergence study for the Cahn-Hilliard equation with no-flux BCs |
| `myodeEX.m` | Explicit part $f_{\text{ex}}(t,y)$ for the scalar ODE test |
| `myodeIM.m` | Implicit operator $G(y)$ for the scalar ODE test (diagonal sparse matrix) |
| `odeIM_CH.m` | Implicit operator $G(t,u)$ for the Cahn-Hilliard equation |
| `BC_periodic.m` | Periodic BCs for the full system $(L, \text{rhs})$ in nonlinear diffusion |
| `BC_periodic_IMEX_L.m` | Periodic BCs applied to system matrix $L$ only (precomputed once per step) |
| `BC_periodic_IMEX_Rhs.m` | Periodic BCs applied to the RHS vector (called each stage) |
| `BC_CH_noflux_IMEX_L.m` | No-flux BCs for the Cahn-Hilliard system matrix $L$ |
| `BC_CH_noflux_IMEX_Rhs.m` | No-flux BCs for the Cahn-Hilliard RHS |
| `BCCH.m` | Combined no-flux BC function for the full Cahn-Hilliard system |
| `filter_spacing.m` | Generates a non-uniform grid with refined spacing near boundaries |

## Citation

If you use this code, please cite the original paper:  
*Semi-implicit-explicit Runge-Kutta Method for Nonlinear Differential Equations.*

```bibtex
@article{ding2025semi,
  title={Semi-implicit-explicit Runge-Kutta method for nonlinear differential equations},
  author={Ding, Lingyun},
  journal={arXiv preprint arXiv:2504.09969},
  year={2025}
}
```

---

## Python Translation

Python translation by **Xi He**.

The MATLAB code has been translated to Python. The folder mapping is:

| MATLAB | Python |
|--------|--------|
| `RKmethod/` | `rk_ode/` |
| `test/` | `rk_ode_test/` |
| `tools/` | `pde_solver/diffmat.py`, `pde_solver/to_latex_table.py` |
| `paper/` | `pde_solver/` |

### `rk_ode/` ← `RKmethod/`

| File | MATLAB equivalent | Purpose |
|------|-------------------|---------|
| `semi_imex.py` | `SemiIMEXRungeKutta.m` | Main semi-IMEX solver |
| `imex.py` | `IMEXRungeKutta2/3.m` | IMEX RK2/RK3 |
| `explicit.py` | `RungeKutta4.m` | RK4 reference solver |
| `step.py` | `RKOneStep.m`, `RKOneStepIMEX2.m` | Butcher tableau and single-step functions |
| `problems.py` | — | Predefined test problems and error metrics |
| `utils.py` | *(no equivalent)* | Python-specific shape and dispatch helpers 

`semiIMEXRungeKutta(..., method=[order, stage, index])` selects a Butcher tableau:

| `method` | Description |
|----------|-------------|
| `[1, 2, 1]` | 1st-order, 2-stage |
| `[2, 2, 1]` | 2nd-order, 2-stage (Table 2) |
| `[2, 3, 4]` | 2nd-order, 3-stage, variant 4 (Table 5) |
| `[2, 3, 3]` | 2nd-order, 3-stage, variant 3 (Table 7) |
| `[3, 4, 1]` | 3rd-order, 4-stage (Table 8) |
| `[3, 5, 1]` | 3rd-order, 5-stage (Table 9) |

#### Note on the revised helper path in `semi_imex.py`

A small helper revision was added in `semi_imex.py` for scalar ODE problems.

What changed:

- a helper `_scalarize_value(x)` was added to safely convert scalar-like outputs such as `()`, `(1,)`, or `(1,1)` into one Python float
- `RKOneStep(...)` now uses a scalar-only fast path when the state size is 1 and no boundary-condition callback is present

Why this was added:

- the scalar convergence benchmark used for the paper is a true 1D ODE
- the generic solver was treating it like a `1x1` matrix, introducing floating-point noise at fine time grids

What the scalar path does:

- stores stage values, explicit terms, and implicit products as plain floats
- solves the scalar implicit stage by direct division instead of calling `np.linalg.solve` on a `1x1` matrix

Vector-valued ODEs and PDE systems still use the generic path.

### `rk_ode_test/` ← `test/`

| File | MATLAB equivalent | Purpose |
|------|-------------------|---------|
| `main.py` | simple test script | Basic sanity check: RK4 on $y'=y$, compares to $e^t$ |
| `test_convergence.py` | convergence test scripts | Convergence rates for semi-IMEX methods; generates LaTeX tables |

### `pde_solver/` ← `tools/` + `paper/`

| File | MATLAB equivalent | Purpose |
|------|-------------------|---------|
| `diffmat.py` | `diffmat.m` | Finite-difference differentiation matrices $[D_0, \ldots, D_m]$ |
| `to_latex_table.py` | `to_latex_convergence_table.m` | Generates publication-ready LaTeX convergence tables |
| `bc.py` | `BC_*.m` files | Boundary condition callbacks for periodic and no-flux BCs |
| `test_scalar_ode_convergence.py` | `paper_scalar_equation.m` | Reproduces **Table 11** |
| `paper_scalar_function.py` | `paper_scalar_equation.m` | Alternative script for the same scalar equation |
| `test_scalar_second_block.txt` | — | Output: errors and rates for 2nd-order methods (Table 11, top block) |
| `test_scalar_third_block.txt` | — | Output: errors and rates for 3rd-order methods (Table 11, bottom block) |
| `test_nonlinear_diffusion.txt` | — | Output: convergence results for the nonlinear diffusion PDE test |

To reproduce Table 11:

```bash
python pde_solver/test_scalar_ode_convergence.py
```

### Setup

```bash
python3 -m venv .venv
.venv/bin/pip install numpy scipy matplotlib
```

### Quick Start

```python
import numpy as np
from rk_ode import semiIMEXRungeKutta

def f_ex(t, y): return np.cos(t) * y      # explicit part
def f_im(t, y): return -y + np.cos(t)     # implicit part (linear in y)

tt = np.linspace(0, 0.5, 1025)            # 1024 steps
y0 = np.array([1.0])

y = semiIMEXRungeKutta(f_ex, f_im, tt, y0, method=[3, 4, 1])  # 3rd-order, Table 8
```
