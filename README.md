# RFIM-Inspired Image Denoising with MCMC

This project explores how **Markov Chain Monte Carlo (MCMC)** methods—specifically **Metropolis-Hastings** and **Gibbs sampling**—can be applied to **binary image denoising**, by modeling pixels as spins in an **Ising-like Markov Random Field (MRF)**.

When the observed image acts as an **external field**, the problem becomes analogous to the **Random Field Ising Model (RFIM)** from statistical physics.  
This connection allows the study of **phase transitions** between ordered and disordered states as image noise increases.  
An additional theoretical appendix extends this analysis using **Replica and Approximate Message Passing (AMP)** methods to characterize information-theoretic phase transitions.

---

## 🔸 Overview

- **Models**
  - *Ising prior:* pairwise interaction enforcing local alignment.
  - *Anisotropic extension:* different couplings along horizontal, vertical, and diagonal directions.
  - *Posterior model:* includes a likelihood term acting as an external field derived from the observed noisy image.

- **Algorithms**
  - **Metropolis sampler:** stochastic spin-flip acceptance based on local energy changes.
  - **Gibbs sampler:** exact conditional resampling at each site.
  - **Simulated annealing:** optional cooling schedule to approximate MAP estimates.

- **Applications**
  - Binary image denoising
  - Study of noise-induced phase transitions
  - Illustration of physical inference concepts (energy minimization ↔ MAP estimation)

---

![Graphical model of the MRF for image denoising](MRF.png)

*Each node \(x_i\) is a latent clean pixel, and \(y_i\) its observed noisy counterpart.
Edges encode spatial interactions (Ising prior), while vertical links represent the external field induced by the data likelihood.*

## 🔸 Phase Transition Interpretation

As the noise level increases, the external field (the observed image) becomes inconsistent, and the system undergoes a qualitative transition:

- **Low noise:** spins align with the field, preserving image structure.  
- **Critical noise:** partial disorder emerges—analogous to the RFIM critical regime.  
- **High noise:** global order collapses; denoising fails.

This mirrors the **order–disorder transition** known from the **Random Field Ising Model**, as discussed in:

> M. Mézard and A. Montanari, *Information, Physics, and Computation*, Oxford University Press, 2009.

---

## 🔸 Repository Structure
rfim-image-denoising

| ├── main.py # Entry point for experiments  
| ├── Gibbs.py # Gibbs sampler for isotropic MRF  
| ├── Metropolis.py # Metropolis-Hastings sampler  
| ├── AGibbs.py # Anisotropic Gibbs sampler  
| ├── AMetropolis.py # Anisotropic Metropolis sampler  
│
├── results

│ ├── Results without external field  
│ ├── Result with Image Prior  
│
├── replica_phase_transitions

│ ├── Exercise1.pdf  
│ ├── Exercise2.pdf  
│ ├── Exercise3.pdf  
│ └── README.md  ← Theoretical appendix: replica & AMP analysis  
│
├── paper.pdf  
│
└── README.md  

---

## 🔸 Theoretical Extension — Replica Phase Transitions

The folder `replica_phase_transitions/` provides a **theoretical complement** to the numerical MCMC experiments.  
It summarizes results derived using **Replica**, **Cavity**, and **Approximate Message Passing (AMP)** techniques, following:

> F. Krzakala and L. Zdeborová, *Statistical Physics Methods in Optimization and Machine Learning*, Collège de France, 2022.

These methods analyze the same kind of **order–disorder transition** observed empirically in the denoising experiments — but from a **mean-field analytical perspective**.  
They predict **information-theoretic phase boundaries** separating recoverable from unrecoverable regimes in inference problems under random noise.

The analogy with the RFIM becomes explicit:  
the external field (the noisy image) biases each site, and as disorder increases, the system crosses a transition where consistent reconstruction becomes impossible.  
This provides a direct theoretical link between **spin-glass physics** and **statistical inference**.

---

## 🔸 Usage

```bash
python main.py
```
The script loads a binary image, adds synthetic noise, and applies both Gibbs and Metropolis samplers for comparison.
Parameters such as temperature, coupling constants, and number of iterations can be tuned directly in main.py.
