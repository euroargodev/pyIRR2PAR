import numpy as np
from pyIRR2PAR.JTan2025_Ed_380_443_490_555 import JTan_PAR_from_Ed
from pyIRR2PAR.JPitarch_Ed_380_443_490_555 import JPitarch_PAR_from_Ed


def IRR2PAR(Ed, z):
    """
    Combines the results of the PAR models from
        -   JTan_2025:  Tan, J., Frouin, R., Leymarie, E., Mitchell, B.G., 2025. Modeling underwater photosynthetically available radiation
        profiles from biogeochemical Argo floats using multi-spectral irradiance measurements. Opt. Express, OE 33, 44355–44377.
        https://doi.org/10.1364/OE.566083

        -   JPitarch_2025 : Pitarch, J., Leymarie, E., Vellucci, V., Massi, L., Claustre, H., Poteau, A., Antoine, D., Organelli, E., 2025.
        Accurate estimation of photosynthetic available radiation from multispectral downwelling irradiance profiles.
        Limnology and Oceanography: Methods n/a. https://doi.org/10.1002/lom3.10673

    Returns two separate lists for each model's results, and a vector of the mean between mean_par and PAR (excluding NA/NaN values).

    Parameters:
    Ed -- 2D array with columns [Ed380, Ed443, Ed490, Ed555] in W/m2/nm
    z  -- depth

    Returns:
        - Two lists: first for JTan_2025 results, second for JPitarch_2025 results,
        - mean_PAR : vector of the mean between mean_par and PAR (excluding NA/NaN values)
        - mean_uncertainty : vector of uncertainties

    """
    # Call JTan_PAR_from_Ed from JTan_2025 with input in uW/Cm2
    JTan_PAR, JTan_e = JTan_PAR_from_Ed(np.column_stack((100 * Ed, z)))
    JTan_PAR_results = [JTan_PAR, JTan_e]

    # Remove negative values for JPitarch_2025
    Ed[Ed < 0] = np.nan

    # Modify null and negative depth values for JPitarch_2025, as it uses a logarithmic transformation of z.
    z[z <= 0] = 1e-6

    # Call JPitarch_PAR_from_Ed from JPitarch_2025
    JPitarch_PAR, JPitarch_PAR_b, JPitarch_ep50, JPitarch_IQR_ep = JPitarch_PAR_from_Ed(
        Ed, z
    )

    # Stack JTan_PAR and JPitarch_PAR into a 2D array
    values = np.vstack((JTan_PAR, JPitarch_PAR))

    # Calculate the mean, ignoring NaN values
    mean_PAR = np.nanmean(values, axis=0)

    # uncertainties
    # Convert IQR_ep to absolute uncertainty
    JPitarch_e = JPitarch_PAR * (JPitarch_IQR_ep / 100)
    JPitarch_PAR_results = [
        JPitarch_PAR,
        JPitarch_PAR_b,
        JPitarch_ep50,
        JPitarch_IQR_ep,
        JPitarch_e,
    ]

    # Calculate mean_uncertainty based on NaN conditions
    mean_uncertainty = np.where(
        np.isnan(JPitarch_PAR),  # If JPitarch_PAR is NaN
        np.where(
            np.isnan(JTan_PAR),  # If JTan_PAR is also NaN
            np.nan,  # mean_uncertainty is NaN
            JTan_e,  # mean_uncertainty is JTan_e
        ),
        np.where(
            np.isnan(JTan_PAR),  # If JTan_PAR is NaN
            JPitarch_e,  # mean_uncertainty is JPitarch_e
            # np.sqrt(JTan_e ** 2 + JPitarch_e ** 2) / 2  # If neither is NaN, use  LPU (Low of Propagation of Uncertainties ) method assuming no correlation
            (JTan_e + JPitarch_e)
            / 2,  # If neither is NaN, use  LPU (Low of Propagation of Uncertainties ) method assuming full correlation
        ),
    )

    return JTan_PAR_results, JPitarch_PAR_results, mean_PAR, mean_uncertainty
