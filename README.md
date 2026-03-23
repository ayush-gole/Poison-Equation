# Poisson Equation — 1D Solver

Solving the 1D Poisson equation using **Gaussian elimination** with partial pivoting, implemented from scratch in Python. Includes a Jupyter notebook with inline outputs.

---

## What it does

Solves ∇²φ = −ρ/ε on a 1D domain [0, 1] with Dirichlet boundary conditions (φ(0) = 0, φ(1) = 1). The equation is discretized using second-order finite differences, which turns it into a tridiagonal linear system Ax = b that's then solved using a hand-written Gaussian elimination.

Two cases are compared:
- **ρ = constant** — uniform charge distribution
- **ρ = 0** — no source term (Laplace equation), which should give a straight line between the boundary values

No `numpy.linalg.solve` or any built-in linear solver is used — the entire Gaussian elimination with partial pivoting is written manually.

---

## How to Run

**Python script:**
```bash
python Poison_Equation.py
```

**Notebook:**
Just open `Poison_Equation.ipynb` in Jupyter — outputs are already there.

---

## Parameters

All near the top of the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 50 | Number of grid points |
| `phi_0` | 0.0 | Boundary value at x = 0 |
| `phi_1` | 1.0 | Boundary value at x = 1 |
| `epsilon` | 1.0 | Permittivity |
| `rho_constant` | 1.0 | Charge density for case 1 |

---

## Tech Stack

- **Python** — `numpy`, `matplotlib`
- **Jupyter** — notebook with inline outputs

---

## Physics Background

The Poisson equation — ∇²φ = −ρ/ε — shows up in electrostatics, fluid mechanics, and heat conduction among other things. In 1D with uniform grid spacing h, the second derivative is approximated as (φᵢ₊₁ − 2φᵢ + φᵢ₋₁)/h², which gives a tridiagonal system with −2 on the diagonal and 1 on the off-diagonals. The boundary conditions modify the first and last rows of the RHS vector. Gaussian elimination with partial pivoting then solves this system directly — no iteration needed.
