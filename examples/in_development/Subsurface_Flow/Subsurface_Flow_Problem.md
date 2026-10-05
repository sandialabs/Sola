# Optimal injection control for pressure management

## Use case

This mimics controlling injection or production during a pressure test, hydraulic characterization experiment, geothermal stimulation, or CO$_2$/water injection pilot. The goal is to create a desired pressure perturbation while avoiding excessive pressure.

---
## Control variable
Distributed or localized injection rate: $q(x)$

---
## High-fidelity model: pressure-dependent permeability

Let $p(x)$ be pressure. Use
$$ -\frac{d}{dx} \left[ \frac{k_{\mathrm{HF}}(p)}{\mu} p' \right] = q(x), $$
with

$$ k_{\mathrm{HF}}(p) = k_0 \exp\left(\alpha(p-p_0)\right). $$

This captures permeability change due to effective stress, fracture aperture change, or compaction.

---
## Low-fidelity model

Use a small-pressure-change approximation:
$$ k_{\mathrm{LF}}(p) = k_0 \left[ 1+\alpha(p-p_0) +\frac{1}{2}\alpha^2(p-p_0)^2 \right]. $$
Then
$$ -\frac{d}{dx} \left[ \frac{k_{\mathrm{LF}}(p)}{\mu}p' \right] = q(x). $$

Both models are nonlinear.

---
## Objective

A simple pressure-targeting optimal control problem is
$$ \min_{q,p} \frac{1}{2} \int_0^1 \left(p(x)-p_d(x)\right)^2\,dx + \frac{\beta}{2} \int_0^1 q(x)^2\,dx, $$
subject to the PDE. Here, $p_d(x)$ is a desired pressure profile.
