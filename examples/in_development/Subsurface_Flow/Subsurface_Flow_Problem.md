# Optimal injection control for pressure management

## Use case

This mimics controlling injection or production during a pressure test, hydraulic characterization experiment, geothermal stimulation, or CO$_2$/water injection pilot. The goal is to create a desired pressure perturbation while avoiding excessive pressure.

The implemented example is intentionally discrepancy-rich. The high-fidelity model includes a localized, pressure-triggered leakoff/storage sink. The low-fidelity model uses the same exponential pressure-dependent permeability but omits this nonlinear sink. This makes the low-fidelity optimum noticeably suboptimal when evaluated with the high-fidelity model.

---
## Control variable
Distributed or localized injection rate: $q(x)$

---
## High-fidelity model: pressure-dependent permeability with nonlinear leakoff

Let $p(x)$ be pressure. Use
$$ -\frac{d}{dx} \left[ \frac{k_{\mathrm{HF}}(p)}{\mu} p' \right]
+ \ell(x)(p-p_0)^3 = q(x), $$
with

$$ k_{\mathrm{HF}}(p) = k_0 \exp\left(\alpha(p-p_0)\right), $$

and

$$ \ell(x)(p-p_0)^3, \qquad
\ell(x)=\ell_0\exp\left[-\left(\frac{x-x_c}{w_c}\right)^2\right]. $$

The cubic sink represents pressure-triggered leakoff or fracture storage that is strongly nonlinear but stabilizing.

---
## Low-fidelity model

Use the same exponential pressure-dependent permeability, but ignore the localized leakoff/storage term:
$$ k_{\mathrm{LF}}(p) = k_0 \exp\left(\alpha(p-p_0)\right). $$
Then
$$ -\frac{d}{dx} \left[ \frac{k_{\mathrm{LF}}(p)}{\mu}p' \right] = q(x). $$

Both models are nonlinear, but the low-fidelity model misses the localized cubic-in-pressure sink. This keeps the permeability physics consistent between fidelities and isolates the discrepancy to the omitted leakoff/storage mechanism.

---
## Objective

A simple pressure-targeting optimal control problem is
$$ \min_{q,p} \frac{1}{2} \int_0^1 \left(p(x)-p_d(x)\right)^2\,dx + \frac{\beta}{2} \int_0^1 q(x)^2\,dx, $$
subject to the PDE. Here, $p_d(x)$ is a desired pressure profile.

In the MATLAB example, the default target is
$$p_d(x) = p_0 + 1.00\sin(\pi x) + 0.25\sin(2\pi x),$$
with default sensitivity $\alpha = 0.5$, localized leakoff coefficient $\ell_0=30$, leakoff center $x_c=0.55$, leakoff width $w_c=0.25$, and regularization $\beta = 10^{-6}$.
