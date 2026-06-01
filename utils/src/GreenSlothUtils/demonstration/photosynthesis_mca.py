import copy

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mxlpy import Model, mca


def create_mca_fig(
    model: Model,
    coeff_psii: str | None,
    coeff_psi: str | None,
    coeff_rubisco: str | None,
    coeff_cytb6f: str | None,
    coeff_atp_synthase: str | None,
    rubp: str | None,
    co2: str | None,
    pq: str | None,
    pc: str | None,
    atp: str | None,
    nadph: str | None,
    v_rubisco: str | None,
    v_psii: str | None,
    v_psi: str | None,
    v_cytb6f: str | None,
    v_atp_synthase: str | None,
) -> tuple[plt.Figure, plt.Axes]:
    """Create curated MCA figure from MxLpy model.

    Create a curated MCA figure from a MxLpy model using the response coefficients of several aspects of photosynthesis. If the model does not contain the given parameters or variables, they will not be plotted and the corresponding row/column will be faded out.

    The chosen MCA will look like this:

            | PSII | PSI | Rubisco |
    RuBP    |
    CO2     |

    Args:
        model (Model): MxLpy model to perform MCA on.
        coeff_psii (str | None): Name of response coefficient for PSII in the MxLpy model.
        coeff_psi (str | None): Name of PSI response coefficient in the MxLpy model.
        coeff_rubisco (str | None): Name of Rubisco response coefficient in the MxLpy model.
        coeff_cytb6f (str | None): Name of Cytb6f response coefficient in the MxLpy model.
        coeff_atp_synthase (str | None): Name of ATP Synthase response coefficient in the MxLpy model.
        rubp (str | None): Name of RuBP representation in the MxLpy model.
        co2 (str | None): Name of CO2 representation in the MxLpy model.
        pq (str | None): Name of PQ representation in the MxLpy model.
        pc (str | None): Name of PC representation in the MxLpy model.
        atp (str | None): Name of ATP representation in the MxLpy model.
        nadph (str | None): Name of NADPH representation in the MxLpy model.
        v_rubisco (str | None): Name of Rubisco flux in the MxLpy model.
        v_psii (str | None): Name of PSII flux in the MxLpy model.
        v_psi (str | None): Name of PSI flux in the MxLpy model.
        v_cytb6f (str | None): Name of Cytb6f flux in the MxLpy model.
        v_atp_synthase (str | None): Name of ATP Synthase flux in the MxLpy model.

    Returns:
        tuple[plt.Figure, plt.Axes]: Figure and Axes of the MCA plot.
    """
    model = copy.deepcopy(model)

    # TODO Find other response coefficients to plot
    # TODO Find other variables and fluxes to plot

    # Create list of parameters to scan if there are not None
    to_scan = [
        i
        for i in [
            coeff_psii,
            coeff_psi,
            coeff_rubisco,
            coeff_cytb6f,
            coeff_atp_synthase,
        ]
        if i is not None
    ]

    # Do MCA of selected response coefficients
    variables, fluxes = mca.response_coefficients(
        model, to_scan=to_scan, disable_tqdm=True
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # Point custom names to appropriate given response coefficients in MCA results
    plot_columns = {
        "PSII": coeff_psii,
        "PSI": coeff_psi,
        "Rubisco": coeff_rubisco,
        "Cytb6f": coeff_cytb6f,
        "ATP Synthase": coeff_atp_synthase,
    }

    # Point custom names to appropriate given variables in MCA results
    plot_vars_index = {
        "RuBP": rubp,
        "CO2": co2,
        "PQ": pq,
        "PC": pc,
        "ATP": atp,
        "NADPH": nadph,
    }

    # Point custom names to appropriate given fluxes in MCA results
    plot_flux_index = {
        r"$v_{\text{PSII}}$": v_psii,
        r"$v_{\text{PSI}}$": v_psi,
        r"$v_{\text{rubisco}}$": v_rubisco,
        r"$v_{\text{Cytb6f}}$": v_cytb6f,
        r"$v_{\text{ATP Synthase}}$": v_atp_synthase,
    }
    
    for point, var_str in plot_vars_index.items():
        if var_str not in model.get_variable_names():
            plot_vars_index[point] = None
    
    plot_vars = variables.loc[[i for i in plot_vars_index.values() if i is not None]].copy()
    plot_vars = plot_vars.rename(
        columns={v: k for k, v in plot_columns.items()},
        index={v: k for k, v in plot_vars_index.items()},
    )
    
    plot_fluxes = fluxes.loc[[i for i in plot_flux_index.values() if i is not None]].copy()
    plot_fluxes = plot_fluxes.rename(
        columns={v: k for k, v in plot_columns.items()},
        index={v: k for k, v in plot_flux_index.items()},
    )

    # Add rows and columns with NaN values for variables not in MCA results
    for i in plot_columns.keys():
        if i not in plot_vars.columns:
            plot_vars[i] = np.nan
        if i not in plot_fluxes.columns:
            plot_fluxes[i] = np.nan

    for i in plot_vars_index.keys():
        if i not in plot_vars.index:
            plot_vars.loc[i, :] = np.nan
    for i in plot_flux_index.keys():
        if i not in plot_fluxes.index:
            plot_fluxes.loc[i, :] = np.nan

    # Specific order of rows and columns
    plot_vars = plot_vars[plot_columns.keys()]
    plot_vars = plot_vars.loc[plot_vars_index.keys()]
    plot_fluxes = plot_fluxes[plot_columns.keys()]
    plot_fluxes = plot_fluxes.loc[plot_flux_index.keys()]
    
    
    # Plot heatmap of MCA results
    cmap_vars = copy.copy(plt.get_cmap("YlGnBu_r"))

    colors = ["#8a0100", "#f4f5da", "#010b7c"] 
    cmap_name = "Vars"
    cmap_vars = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)
    cmap_vars.set_bad(color="white")
    max_coeff = np.nanmax(np.abs(plot_vars.values))
    im1 = ax1.imshow(plot_vars.values, cmap=cmap_vars, interpolation="nearest", vmin=-max_coeff, vmax=max_coeff)
    fig.colorbar(im1, ax=ax1)
    ax1.set_title("Variables")
    ax1.set_xticks(np.arange(len(plot_vars.columns)), labels=plot_vars.columns, rotation=45, ha="right")
    ax1.set_xticklabels(plot_vars.columns)
    ax1.set_yticks(np.arange(len(plot_vars.index)))
    ax1.set_yticklabels(plot_vars.index)
    
    colors = ["#518568", "#f4f5da", "#aa5db8"   ] 
    cmap_name = "Fluxes"
    cmap_fluxes = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)
    cmap_fluxes.set_bad(color="white")
    max_coeff = np.nanmax(np.abs(plot_fluxes.values))
    im2 = ax2.imshow(plot_fluxes.values, cmap=cmap_fluxes, interpolation="nearest", vmin=-max_coeff, vmax=max_coeff)
    fig.colorbar(im2, ax=ax2)
    ax2.set_title("Fluxes")
    ax2.set_xticks(np.arange(len(plot_fluxes.columns)), labels=plot_fluxes.columns, rotation=45, ha="right")
    ax2.set_xticklabels(plot_fluxes.columns)
    ax2.set_yticks(np.arange(len(plot_fluxes.index)))
    ax2.set_yticklabels(plot_fluxes.index)

    # Set axis labels and if values are NaN set alpha of text to 0.3
    for text in ax1.get_yticklabels():
        if plot_vars_index[text.get_text()] is None:
            text.set_color("#b3b3b3ff")
            
    for text in ax2.get_yticklabels():
        if plot_flux_index[text.get_text()] is None:
            text.set_color("#b3b3b3ff")
    
    for ax in [ax1, ax2]:
        for text in ax.get_xticklabels():
            if plot_columns[text.get_text()] is None:
                text.set_color("#b3b3b3ff")
                
    plt.tight_layout()

    return fig, (ax1, ax2)