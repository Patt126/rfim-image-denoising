import numpy as np
import random


def Metropolis(Sample, height, width, ITERA, alpha, Beta, T,
               Yobs=None, lam=0.0, anneal=False, T0=None, Tf=None):
    """
    Metropolis-Hastings sampler for a binary MRF / Ising model with optional simulated annealing.

    Energy convention:
        H(x) = α * Σ_i x_i
             + β * Σ_<i,j> x_i x_j
             - λ * Σ_i y_i x_i
    and configurations are weighted by exp( -H / T ).

    Interpretation:
    - β < 0  → alignment / smoothing (ferromagnetic).
    - β > 0  → anti-alignment (checkerboard).
    - λ > 0  → agreement with observed data y_i.
    - T      → temperature; controls acceptance of uphill moves.

    Arguments
    ---------
    Sample : 1D numpy array
        Current configuration (±1 spins) including a 1-pixel border.
        Updated in place.

    height, width : int
        True image dimensions excluding border.

    ITERA : int
        Number of full sweeps.

    alpha : float
        External field weight.

    Beta : float
        Coupling constant between neighbours.

    T : float
        Base temperature (used if anneal=False).

    Yobs : 1D numpy array, optional
        Observed noisy image in {-1,+1}. If None → prior-only mode.

    lam : float
        Data fidelity strength λ. Controls adherence to observed data.

    anneal : bool
        If True, applies exponential temperature decay between T0 and Tf.

    T0 : float, optional
        Starting temperature for annealing. Defaults to current T.

    Tf : float, optional
        Final temperature for annealing. Defaults to 0.1 * T0.

    Returns
    -------
    Out_inner : 2D numpy array (height, width)
        Final spins cropped to remove border.
    """

    H = height + 2
    W = width + 2
    Nsites = height * width

    # Fallback: if no observed image, reuse Sample (lam=0 → pure prior)
    if Yobs is None:
        Yobs = Sample

    # Annealing setup
    if anneal:
        if T0 is None:
            T0 = T
        if Tf is None:
            Tf = 0.1 * T0

    # Interior indices (skip padded border)
    indices = np.zeros(Nsites, dtype=int)
    t = 0
    for i in range(1, H - 1):
        base = i * W
        for j in range(1, W - 1):
            indices[t] = base + j
            t += 1

    # Main Metropolis loop
    for sweep in range(ITERA):

        # Simulated Annealing Schedule 
        if anneal:
            # Exponential decay of temperature across sweeps
            T_eff = T0 * ((Tf / T0) ** (sweep / max(1, ITERA - 1)))
        else:
            T_eff = T

        np.random.shuffle(indices)

        for k in indices:
            s = Sample[k]    # current spin (±1)
            yk = Yobs[k]     # observed pixel (±1)

            # 4-neighbour contribution
            nb_sum = (
                Sample[k - 1] +
                Sample[k + 1] +
                Sample[k - W] +
                Sample[k + W]
            )


            E_cur = alpha * s + Beta * s * nb_sum - lam * yk * s
            s_new = -s
            E_new = alpha * s_new + Beta * s_new * nb_sum - lam * yk * s_new

            dE = E_new - E_cur

            # Metropolis acceptance rule
            if dE <= 0:
                Sample[k] = s_new
            else:
                if random.random() < np.exp(-dE / T_eff):
                    Sample[k] = s_new

    # Convert back to 2D and remove padding
    Out_full = np.zeros((H, W))
    for i in range(1, H - 1):
        base = i * W
        for j in range(1, W - 1):
            Out_full[i, j] = Sample[base + j]

    Out_inner = Out_full[1:-1, 1:-1]
    return Out_inner
