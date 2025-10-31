# Replica Phase Transitions — Theoretical Appendix

This folder extends the MCMC image denoising project with a theoretical perspective from the **statistical physics of inference**, linking spin-based image models to the analysis of phase transitions in high-dimensional recovery problems.

---

## Overview

The material builds on the framework introduced in  
*Florent Krzakala* and *Lenka Zdeborová*,  
_Statistical Physics Methods in Optimization and Machine Learning_ (Collège de France, 2022).

It explores the **Replica method** and **Approximate Message Passing (AMP)** algorithms used to study disordered systems and inference problems.  
In these models, the interplay between noise, coupling strength, and prior structure gives rise to **information-theoretic phase transitions**—sharp thresholds separating recoverable, partially recoverable, and unrecoverable regimes.

---

## Main Concepts

- **Replica Method:** Analytical tool to compute the typical free energy and mean-square error of inference problems under random disorder.  
- **Cavity and AMP Approaches:** Algorithmic counterparts of replica predictions, enabling efficient iterative inference close to theoretical limits.  
- **Phase Transitions:** Critical thresholds at which signal reconstruction abruptly fails as noise or system size crosses a boundary.  
- **Random Field Ising Analogy:** The denoising MRF model behaves as a spatial Random Field Ising Model, where the noisy image acts as an external field breaking spin symmetry.

---

## Connection to the Denoising Project

In the MCMC denoising experiments, increasing noise intensity produces a transition from ordered (recoverable) configurations to disordered (unrecoverable) ones.  
This empirical observation parallels the **theoretical phase boundaries** derived from replica symmetry in mean-field models.  

Both frameworks describe the same phenomenon:
> A competition between interaction strength (coupling) and disorder (noise) that determines the limit of reliable inference.

The included notes illustrate how analytical tools like the Replica trick predict such thresholds, providing a theoretical foundation for the observed dynamics in Ising-type denoising models.

---

## Included Files

- `Exercise1.pdf` — Replica-symmetric derivation for a binary inference model  
- `Exercise2.pdf` — Derivation and fixed-point analysis of Approximate Message Passing (AMP)  
- `Exercise3.pdf` — Phase diagram and recoverability thresholds  

---

## Reference
F. Krzakala and L. Zdeborová,  
*Statistical Physics Methods in Optimization and Machine Learning*,  
Collège de France, 2022.
