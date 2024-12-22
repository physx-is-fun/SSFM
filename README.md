# Laser pulse propagation in optical fiber

This notebook shows a simple, scalar implementation of the [Split-Step Fourier Method](https://en.wikipedia.org/wiki/Split-step_method) for solving the [Nonlinear Schrödinger Equation](https://en.wikipedia.org/wiki/Nonlinear_Schrödinger_equation).

$\frac{\partial A}{\partial z}=-\frac{\alpha}{2}A +i \frac{\beta_2}{2} \frac{\partial^2 A}{\partial t^2} + \frac{\beta_3}{6} \frac{\partial^3 A}{\partial t^3} -i \gamma[(|A|^2A) -\frac{i}{\omega_0}\frac{\partial}{\partial t}(|A|^2A) -T_R A\frac{\partial}{\partial t}(|A|^2)]$

This nonlinear partial differential equation models how the envelope and phase of light pulse changes when propagating through a single mode optical fiber, when taking power attenuation ($\alpha$), group velocity dispersion ($\beta_2$), third order dispersion ($\beta_3$), waveguide nonlinearity ($\gamma$) causing self-phase modulation (SPM), self-steepening (SS) and Stimulated Raman-response (SRS) ($T_R$) into account. Note that the effects of SS and and SRS are not included (yet) in the code.

## References

    Ole Krarup. (2024). NLSE-vector-solver. Retrieved from https://github.com/OleKrarup123/NLSE-vector-solver
    Ole Krarup. (2022). The Nonlinear Schrödinger Equation solved in python!. Retrieved from https://colab.research.google.com/drive/1XyLYqrohf5GL6iFSrS6VlHoj_eSm-kAG?usp=sharing