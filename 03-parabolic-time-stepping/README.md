Implementation of a solver for the non-stationary parabolic diffusion equation using a space semi-discrete Galerkin approach coupled with the $\theta$-method family for time integration. The matrix differential system $P\dot{u} + Au = b$ incorporates the analytically computed mass matrix $P$ and the stiffness matrix $A$. The project compares Backward Euler ($\theta=1$), Crank-Nicolson ($\theta=0.5$), and Forward Euler ($\theta=0$), assessing unconditional stability properties against the stringent parabolic CFL constraint ($\Delta t \le C h^2$) inherent to the explicit scheme. Finally, the analysis quantifies the global space-time error measured in the $L^2(0,T; L^2(\Omega))$ norm, validating the optimal temporal balancing $\Delta t = \mathcal{O}(h)$ for second-order schemes and $\Delta t = \mathcal{O}(h^2)$ for first-order schemes.

[complete report](./Numerical_Methods_for_PDEs_03.pdf/)

Codes --> 'Project3' folder
