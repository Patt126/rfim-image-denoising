from math import exp
import random
import numpy as np


def AGibbs(Sample, height, width, ITERA,
           alpha, Beta_h, Beta_v, Beta_d1, Beta_d2,
           T, Yobs=None, lam=0.0, anneal=False, T0=None, Tf=None):
    """
    Anisotropic Gibbs sampler for binary image denoising via an Ising MRF posterior,
    extended to separate horizontal, vertical, and diagonal couplings,
    and with optional simulated annealing and likelihood term.

    Energy convention:
        E(X|Y) = Σ_k [
            α * X_k
            + (β_h * X_k * Σ_horiz N_h)
            + (β_v * X_k * Σ_vert N_v)
            + (β_d1 * X_k * Σ_diag1 N_d1)
            + (β_d2 * X_k * Σ_diag2 N_d2)
            - λ * Y_k * X_k
        ]
    and states are weighted by exp(-E / T).

    Notes:
    - β_h, β_v, β_d1, β_d2 < 0 → alignment in corresponding directions.
    - λ > 0 → attraction to observed data (data fidelity).
    - T controls randomness; can decay exponentially if anneal=True.

    Parameters
    ----------
    Sample : 1D numpy array
        Current spin configuration in {-1,+1} including a 1-pixel border.

    height, width : int
        True image size excluding the 1-pixel border.

    ITERA : int
        Number of full sweeps.

    alpha : float
        External field parameter.

    Beta_h, Beta_v, Beta_d1, Beta_d2 : float
        Coupling constants for horizontal, vertical, and diagonal neighbours.

    T : float
        Base temperature (used if anneal=False).

    Yobs : 1D numpy array, optional
        Observed noisy image (±1). If None, defaults to prior-only mode.

    lam : float
        Data fidelity term λ.

    anneal : bool
        Whether to apply simulated annealing.

    T0, Tf : float, optional
        Starting and final temperatures for annealing. Defaults: T0=T, Tf=0.1*T0.

    Returns
    -------
    SampleOut : 2D numpy array (height x width)
        Final denoised configuration (without padding).
    """

    # original image dimensions
    H = height + 2
    W = width + 2
    Nsites = height * width

    # fallback for missing Yobs
    if Yobs is None:
        Yobs = Sample

    # annealing parameters
    if anneal:
        if T0 is None:
            T0 = T
        if Tf is None:
            Tf = 0.1 * T0

    # precompute interior indices
    indices = np.zeros(Nsites, dtype=int)
    t = 0
    for i in range(1, H - 1):
        for j in range(1, W - 1):
            indices[t] = j + i * W
            t += 1

    # main Gibbs loop
    for sweep in range(ITERA):

        # simulated annealing schedule
        if anneal:
            T_eff = T0 * ((Tf / T0) ** (sweep / max(1, ITERA - 1)))
        else:
            T_eff = T

        random.shuffle(indices)

        for k in indices:
            # observed pixel
            yk = Yobs[k]

            # compute local energy contribution from anisotropic neighbours
            # horizontal: left/right
            # vertical: up/down
            # diagonal1: top-left / bottom-right
            # diagonal2: top-right / bottom-left
            energy = (
                alpha
                + Beta_h * (Sample[k - 1] + Sample[k + 1])
                + Beta_v * (Sample[k - W] + Sample[k + W])
                + Beta_d1 * (Sample[k - W - 1] + Sample[k + W + 1])
                + Beta_d2 * (Sample[k - W + 1] + Sample[k + W - 1])
                - lam * yk  # data fidelity term
            )

            # conditional probability 
            p_plus = 1.0 / (1.0 + exp(2.0 * energy / T_eff))

            # sample from Bernoulli
            Sample[k] = 1 if random.random() < p_plus else -1

    # reshape 1D sample into 2D output without border
    SampleOut = np.zeros((H, W))
    for i in range(1, H - 1):
        for j in range(1, W - 1):
            SampleOut[i][j] = Sample[j + i * W]

    return SampleOut[1:-1, 1:-1]
