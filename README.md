# Solving Burgers' Equation using Spectral Method

## Introduction

In this project, we aim to solve the viscous Burgers' equation using a spectral method. The significance of the viscous Burgers' equation lies in its resemblance to the momentum equation of the Navier-Stokes (NS) equation, albeit without the complexities of three-dimensional flow and pressure-velocity coupling issues.

Turbulence exhibits spectral characteristics, making spectral decomposition an advantageous method for solving turbulence-related problems. A spectral solver decomposes a Partial Differential Equation (PDE) into an infinite sum of Ordinary Differential Equations (ODE) in Fourier space. Each ODE is then solved independently, and the overall solution of the PDE is obtained by summing the solutions in Fourier space.

Given that we truncate the Fourier modes up to a certain wave number ($K_N$), a turbulence model is necessary when the cutoff wave number does not cover the majority of the energy spectrum. Consequently, we introduce a spectral eddy-viscosity model to address this issue when solving the Burgers' equation.

In the appendix, we provide a fast and computationally efficient MATLAB code that avoids using `if/else` conditions.

## Governing Equations

The viscous Burgers' equation without a forcing term is given by:

\[ \partial_t u + u \partial_x u = \frac{1}{Re} \partial_{x x} u \]

The transformation of this equation to Fourier space is already discussed in the class and lecture notes, yielding the final equation:

\[ \partial_t \hat{u}_k + \sum_{k=p+q} \hat{u}_p i q \hat{u}_q = -\nu k^2 \hat{u}_k \quad k=0, \ldots, N \]

where $\nu=\frac{1}{Re}$.

When employing a Large Eddy Simulation (LES) with a spectral Subgrid-Scale (SGS) model, the velocity $\hat{u}$ becomes filtered ($\hat{\bar{u}}$), and the diffusive term is modified as:

\[ \partial_t \hat{\bar{u}}_k + \sum_{k=p+q} \hat{\bar{u}}_p i q \hat{\bar{u}}_q = -k^2 \nu_{eff}(k) \hat{\bar{u}}_k \quad k=0, \ldots, N \]

where $\nu_{\text{eff}}(k) = \nu + \nu_t(k)$, and $\nu_t$ is the turbulent viscosity.

The turbulent viscosity $\nu_t$ is calculated using the eddy-viscosity model proposed by Metais and Lesieur:

\[ 
\left\{
\begin{split}
\nu_t &= \nu_t^{+\infty} \left( \frac{E_{k_N}}{k_N} \right)^{1 / 2} \nu_t^* \left( \frac{k}{k_N} \right) \\
\nu_t^{+\infty} &= 0.31 \frac{5-m}{m+1} \sqrt{3-m} C_K^{-3 / 2} \\
\nu_t^* &= 1+34.5 e^{-3.03(k_N / k)} \\
\end{split}
\right.
\]

where $k_N$ is the cutoff wave number, $E_{k_N}$ is its corresponding energy, $m$ is the slope of the inertial region of the energy spectrum, and $C_K$ is the Kolmogorov constant. For the 1D Burgers equation $m = 2$ and $ C_k \approx 0.4523$.

## Numerical Method

The forward Euler formula can easily discretize equations \ref{eq:1.2} and \ref{eq:1.3}. However, the triadic convection term should be split into a summation in terms of $p$ and $q$ for the sake of coding. For any $0 \leq k \leq N $, $\sum_{k=p+q} = \left.\sum_{p = k-N}^{N}\right|_{q =k-p}$. Therefore, \autoref{eq:1.3} can be discretized as:

\[ \hat{u}_k^{n+1} = \hat{u}_k^{n} + \Delta t \cdot (-\left.\sum_{p = k-N}^{N}\right|_{q =k-p}\hat{u}_p^{n} i q \hat{u}_q^{n} -k^2\nu_{eff}(k) \hat{u}_k^{n}) \quad k=0, \ldots, N \]

At each time step, we use the conjugate property to calculate $\hat{u}_k^{n}$ for $-N \leq k \leq -1$:

\[ \hat{u}_k^{n} = - \hat{u}_{-k}^{n}, \quad k = -N, \ldots, -1 \]

The energy at each time step is calculated using:

\[ \hat{E}_k^{n} = \hat{u}_k^{n} \overline{\hat{u}_k^{n}} \]

## Results and Discussion

Figure~\ref{fig:Ek} presents the energy spectrum for different scenarios. For $Re = 40$, a cutoff wave number of 100 captures the energy spectrum up to the dissipation region. Conversely, the $N = 20$ case is limited to the inertial region, necessitating a turbulent viscosity model.

The two LES cases illustrate the effect of $C_K$ on the spectrum. A too-small value for $C_K$ results in the dissipative region starting from smaller wave numbers, making the flow excessively dissipative, as observed in the $C_K = 0.05$ case.

![Energy spectrum of the steady-state solution of the Burgers equation for different cases.](figs/Ek.eps)

Figure~\ref{fig:para} demonstrates a parametric study on the effect of $Re$ number on the energy spectrum. Increasing $Re$ makes the dissipative region appear at higher wave numbers, while the inertial region remains unaffected. This observation aligns with the expected behavior: for a given $k$, an increase in $Re$ leads to a decrease in viscosity, resulting in reduced diffusion.

For a given $k$: 
\[ Re \uparrow \quad \rightarrow \quad \nu \downarrow \quad \rightarrow \quad \nu \frac{\partial^2 u}{\partial x^2} \downarrow \]

![Effect of $Re$ number on the energy spectrum. Blue curve: $Re = 200$, Orange curve: $Re = 500$.](figs/parametric.jpg)

## Conclusion

In this study, we explored the numerical solution of the viscous Burgers' equation using a spectral element method. The significance of this equation lies in its similarity to the momentum equation of the Navier-Stokes equation, making it a valuable model for understanding turbulence dynamics without the complexities of three-dimensional flow and pressure-velocity coupling.

We found that employing a spectral method offers several advantages for solving turbulence problems, given the spectral nature of turbulence. By decomposing the partial differential equation into an infinite sum of ordinary differential equations in Fourier space, we can efficiently solve for the system's dynamics.

Furthermore, we introduced a spectral eddy-viscosity model to account for unresolved turbulent motions, demonstrating the importance of incorporating turbulence modeling in simulations where the cutoff wave number does not capture the entire energy spectrum.

Our numerical results revealed insights into the behavior of the energy spectrum under different scenarios. Notably, we observed the influence of the Reynolds number on the appearance of the dissipative region in the spectrum, highlighting the importance of Reynolds number considerations in turbulence studies.

In conclusion, our study sheds light on the numerical solution of Burgers' equation and its implications for understanding turbulence dynamics. Future research could explore more sophisticated turbulence models and investigate the applicability of spectral methods to more complex flow configurations.
