# Dynamical Systems
In this research focus is directed upon dynamic robotic perception domains in which communication difficulties require autonomous control. Two of such domains are autonomous underwater vehicles (AUVs) and spacecraft.

## Autonomous Underwater Vehicle
For sake of simplicity to demonstrate control capabilities, it suffices to model an autonomous underwater vehicle (AUV) as a submerged sphere, subject to gravity, buoyancy, fluid drag, and omnidirectional thrust.

### State Equations of Motion
$$
\vec{s} =
\\left[\\begin{matrix}
\dot{\vec{p}} \\\\
\dot{\vec{v}} \\\\
\\end{matrix}\\right]
=
\\left[\\begin{matrix}
\vec{v} \\\\
m\vec{g} - \frac{1}{2} \rho v \vec{v} C_D - \rho V \vec{g} + \frac{T_{max}}{m} \vec{u} \\\\
\\end{matrix}\\right]
$$

### Cost Functional
$$
J (t_0, t_1, \vec{u}(t)) = \alpha \int_{t_0}^{t_1} u(t) dt + (1 - \alpha) \int_{t_0}^{t_1} u^2(t)dt
$$

### Hamiltonian
$$
\mathcal{H} =
\vec{\lambda} \cdot \vec{s} + \alpha u + (1 - \alpha) u^2
$$

### Optimal Thrust Direction
$$
\hat{u}^\star = - \frac{\vec{\lambda}_v}{\lambda_v}
$$

### Optimal Thrust Magnitude
$$
u^\star = \frac{\alpha m - T_{max} \lambda_v}{2 m (\alpha - 1)}
$$

Noting $u \in [0, 1]$ and the singularity in $u^\star$ when $\alpha = 1$, the optimal thrust when $\alpha < 1 $ is
$$
u^\star =
\text{min}\left(
  \text{max}\left(
    \frac{\alpha m - T_{max} \lambda_v}{2 m (\alpha - 1)}, 0
    \right), 1\right)
$$
and when $\alpha = 1$
$$
u^\star =
\begin{cases}
1 & \text{if} & s > 1 \\\\
0 & \text{if} & s < 0 \\\\
\end{cases}
$$
where the "switching function" is
$$
s = \alpha m - T_{max} \lambda_v
$$
This discontinuous control policy with $\alpha = 0$ describes "bang-bang" optimal control, and in this system describes the minimisation of power usage.
