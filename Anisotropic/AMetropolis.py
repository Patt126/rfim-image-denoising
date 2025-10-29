from math import exp
import random
import numpy as np


def AMetropolis(Sample,height,width,ITERA,Alpha,Beta_h, Beta_v, Beta_d1, Beta_d2,T,Yobs=None,lam=0.0,anneal=False,T0=None,Tf=None):
    """
    Anisotropic Metropolis-Hastings sampler for binary MRF denoising.

    State convention:
    - Sample[k] ∈ {-1,+1} is the current latent (clean) pixel value.
    - Yobs[k]   ∈ {-1,+1} is the observed noisy pixel (fixed, does not change).
    - The lattice is stored with a 1-pixel border around the true image.

    Energy convention (anisotropic):
        E(X|Y) = Σ_k [
            Alpha        * X_k
          + Beta_h  * X_k * (left+right)
          + Beta_v  * X_k * (up+down)
          + Beta_d1 * X_k * (diag ↘ / ↖)
          + Beta_d2 * X_k * (diag ↗ / ↙)
          - lam     * Y_k * X_k
        ]

    and states are weighted by exp( -E / T ).

    Notes:
    - Negative Beta_h, Beta_v, Beta_d1, Beta_d2 encourage alignment along their respective directions.
    - lam > 0 encourages agreement with the observed data Yobs.
    - T controls how often we accept energetically worse flips.
    - anneal=True activates simulated annealing: T decays from T0 to Tf.

    Parameters
    ----------
    Sample : 1D numpy array
        Spin configuration (±1) including border. Will be modified in place.

    height, width : int
        True image dimensions (excluding the padding border).

    ITERA : int
        Number of full sweeps.

    Alpha : float
        Single-site bias term.

    Beta_h, Beta_v, Beta_d1, Beta_d2 : float
        Directional couplings:
        - Beta_h  : horizontal (left/right)
        - Beta_v  : vertical (up/down)
        - Beta_d1 : main diagonal (top-left / bottom-right)
        - Beta_d2 : other diagonal (top-right / bottom-left)

    T : float
        Base temperature for Metropolis acceptance if anneal=False.

    Yobs : 1D numpy array, optional
        Observed noisy image in {-1,+1}, same shape/indexing as Sample.
        If None → prior-only version (no data term).

    lam : float
        Data fidelity weight λ. Higher λ forces agreement with Yobs.

    anneal : bool
        If True, enable exponential temperature decay across sweeps.

    T0, Tf : float, optional
        Start and end temperatures for annealing. Defaults: T0=T, Tf=0.1*T0.

    Returns
    -------
    SampleOut_no_border : 2D numpy array of shape (height, width)
        Final spin configuration cropped to remove the border.
    """

    H = height + 2   # total rows including 1-pixel border
    W = width + 2    # total cols including 1-pixel border
    Nsites = height * width  # number of interior pixels

    # If Yobs is not provided, fall back to prior-only behaviour
    if Yobs is None:
        Yobs = Sample

    # Annealing setup
    if anneal:
        if T0 is None:
            T0 = T
        if Tf is None:
            Tf = 0.1 * T0

    # Build list of interior indices (skip padded frame)
    indices = np.zeros(Nsites, dtype=int)
    t = 0
    for i in range(1, H - 1):
        for j in range(1, W - 1):
            indices[t] = j + i * W
            t += 1

    # Main Metropolis loop over sweeps
    for sweep in range(ITERA):

        # Select current effective temperature for this sweep
        if anneal:
            # exponential decay from T0 down to Tf
            T_eff = T0 * ((Tf / T0) ** (sweep / max(1, ITERA - 1)))
        else:
            T_eff = T

        random.shuffle(indices)

        for k in indices:
            # Current spin at site k
            current = Sample[k]  # ±1

            # Proposed flip
            candidate = -current

            # Observed pixel value (data fidelity term)
            yk = Yobs[k]  # ±1

            # Compute local energy contribution before the flip:
            # Minus lam * yk * spin is the likelihood attachment.
            E_old = (
                Alpha * current
                + Beta_h  * (current * Sample[k - 1]       + current * Sample[k + 1])
                + Beta_v  * (current * Sample[k - W]       + current * Sample[k + W])
                + Beta_d1 * (current * Sample[k - W - 1]   + current * Sample[k + W + 1])
                + Beta_d2 * (current * Sample[k - W + 1]   + current * Sample[k + W - 1])
                - lam     * yk * current
            )

            # Compute local energy after flipping that one spin:
            E_new = (
                Alpha * candidate
                + Beta_h  * (candidate * Sample[k - 1]       + candidate * Sample[k + 1])
                + Beta_v  * (candidate * Sample[k - W]       + candidate * Sample[k + W])
                + Beta_d1 * (candidate * Sample[k - W - 1]   + candidate * Sample[k + W + 1])
                + Beta_d2 * (candidate * Sample[k - W + 1]   + candidate * Sample[k + W - 1])
                - lam     * yk * candidate
            )

            Delta_E = E_new - E_old  # energy change for this flip

            # Metropolis acceptance rule with temperature T_eff:
            # If energy goes down (Delta_E < 0) accept immediately.
            # Otherwise accept with prob exp(-Delta_E / T_eff).
            if Delta_E < 0:
                Sample[k] = candidate
            else:
                if random.random() < exp(-Delta_E / T_eff):
                    Sample[k] = candidate

    # Convert flattened lattice with border back to 2D and drop the artificial border
    SampleOut = np.zeros((H, W))
    for i in range(1, H - 1):
        for j in range(1, W - 1):
            SampleOut[i][j] = Sample[j + i * W]

    SampleOut_no_border = SampleOut[1:-1, 1:-1]
    return SampleOut_no_border
