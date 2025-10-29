import numpy as np
import random


def Gibbs(Sample, height, width, ITERA, alpha, Beta, T,
          Yobs=None, lam=0.0, anneal=False, T0=None, Tf=None):
    """
    Gibbs sampler for binary image denoising via an Ising MRF posterior,
    with optional simulated annealing.

    Energy convention:
        E(X|Y) = Σ_k [ α * X_k + β * X_k * Σ_{n∈N(k)} X_n - λ * Y_k * X_k ]
    and configurations are weighted by exp( -E / T ).

    Notes:
    - β < 0 encourages alignment (ferromagnetic smoothing).
    - β > 0 encourages anti-alignment (checkerboard).
    - λ > 0 enforces similarity with observed data Yobs.
    - Temperature T controls randomness; lower → more deterministic.

    Annealing:
    - If anneal=True, temperature decays exponentially from T0 to Tf.

    Arguments
    ---------
    Sample : 1D numpy array
        Current spin configuration (±1) including 1-pixel border.
        Updated in place.

    height, width : int
        True image dimensions (without border).

    ITERA : int
        Number of full sweeps.

    alpha : float
        External field term.

    Beta : float
        Coupling constan.

    T : float
        Base temperature if anneal=False.

    Yobs : 1D numpy array, optional
        Observed noisy image in {-1,+1}.
        If None, prior-only sampling (lam=0).

    lam : float
        Data fidelity parameter (strength of attraction to Yobs).

    anneal : bool
        Enable simulated annealing (temperature decay).

    T0 : float, optional
        Starting temperature for annealing (default = T).

    Tf : float, optional
        Final temperature for annealing (default = 0.1 * T0).

    Returns
    -------
    Out_inner : 2D numpy array
        Final configuration (no padding).
    """

    H = height + 2
    W = width + 2
    Nsites = height * width

    # Fallback for missing observation (pure prior mode)
    if Yobs is None:
        Yobs = Sample

    # Setup annealing parameters
    if anneal:
        if T0 is None:
            T0 = T
        if Tf is None:
            Tf = 0.1 * T0

    # Precompute linear indices for interior pixels (exclude padded border)
    indices = np.zeros(Nsites, dtype=int)
    t = 0
    for i in range(1, H - 1):
        base = i * W
        for j in range(1, W - 1):
            indices[t] = base + j
            t += 1

    # Main Gibbs sampling loop
    for sweep in range(ITERA):

        #Simulated Annealing Schedule
        if anneal:
            # Exponential decay from T0 to Tf
            T_eff = T0 * ((Tf / T0) ** (sweep / max(1, ITERA - 1)))
        else:
            T_eff = T
        np.random.shuffle(indices)

        for k in indices:
            yk = Yobs[k]

            # 4-neighbour sum (standard Ising coupling)
            nb_sum = (
                Sample[k - 1] +
                Sample[k + 1] +
                Sample[k - W] +
                Sample[k + W]
            )

            # α + β * nb_sum - λ * yk
            # Negative β favours alignment.
            h_loc = alpha + Beta * nb_sum - lam * yk

            # Conditional probability for X_k = +1.
            #   p_plus = 1 / (1 + exp( 2*h_loc / T_eff ))
            p_plus = 1.0 / (1.0 + np.exp(2.0 * h_loc / T_eff))

            # Draw new spin according to conditional probability.
            Sample[k] = 1 if random.random() < p_plus else -1

    # Reconstruct 2D array and crop border.
    Out_full = np.zeros((H, W))
    for i in range(1, H - 1):
        base = i * W
        for j in range(1, W - 1):
            Out_full[i, j] = Sample[base + j]

    Out_inner = Out_full[1:-1, 1:-1]
    return Out_inner
