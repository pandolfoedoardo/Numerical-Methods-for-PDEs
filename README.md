# Numerical Methods for Partial Differential Equations

A collection of numerical implementations and convergence studies for elliptic, parabolic, and convection-dominated PDEs using Finite Element Methods (FEM) and finite difference time-stepping.

## Projects Overview

* [**01. Elliptic Problems & FEM Discretization**](./01-elliptic-problems-fem/)
  Weak formulation, stiffness matrix assembly, and empirical $L^2$ / $H^1$ error convergence rates.
* [**02. Parabolic PDEs & Time-Stepping Schemes**](./02-parabolic-time-stepping/)
  Space discretization via Galerkin FEM coupled with implicit/explicit time integrators (Euler, Crank-Nicolson).
* [**03. Convection-Diffusion & Stabilization Techniques**](./03-convection-diffusion-supg/)
  Analysis of Péclet number instabilities and implementation of SUPG (Streamline Upwind Petrov-Galerkin) stabilization.
