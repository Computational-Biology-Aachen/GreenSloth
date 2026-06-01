import copy
from datetime import datetime

import matplotlib.pyplot as plt
import neonutilities as nu
import pandas as pd
from .utils import custom_latex
from mxlpy import Model, Simulator, make_protocol
import matplotlib.dates as mdates
from IPython.display import clear_output
from pysunnoaa import noaa

def strike(text: str) -> str:
    return "".join([f"\u0336{c}" for c in text])

def simulate_again(s: Simulator, day_prtc, time_points) -> None:
    time_points += 100
    if time_points > 1e4:
        return None
    s.simulate_protocol(day_prtc, time_points_per_step=time_points).get_result().default(simulate_again(s, day_prtc, time_points))
    return 
    

def create_day_simulation_fig(
    model: Model,
    pfd: str,
    vc: str | None = None,
    atp: str | None = None,
    nadph: str | None = None,
    fluorescence: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Create a day simulation figure.

    Create a day simulation figure using a MxLpy model and PAR data from NEON (https://www.neonscience.org/). The figure consists of the PAR data as a filled area plot and the model results as line plots. The model results can include the Rubisco carboxylation rate, the ATP/NADPH ratio, and the fluorescence. If the model does not contain the variables, they will not be plotted and "n.a." will be written as the y-axis label.

    Args:
        model (Model): MxLpy model to simulate day with
        pfd (str): Name of the PPFD parameter in the MxLpy model
        vc (str | None, optional): Name of the Rubisco carboxylase activity rate in the MxLpy model. Defaults to None.
        atp (str | None, optional): Name of the ATP variable in the MxLpy model. Defaults to None.
        nadph (str | None, optional): Name of the NADPH variable in the MxLpy model. Defaults to None.
        fluorescence (str | None, optional): Name of the fluorescence variable in the MxLpy model. Defaults to None.

    Returns:
        tuple[plt.Figure, plt.Axes]: Figure and axis of the day simulation plot
    """
    model = copy.deepcopy(model)
    
    # TODO: Check for case if result = None
    # TODO: Seperate Simulation results to check if given name is in variables or fluxes or parameters or readouts or surrogates

    # Load PAR data from NEON at the KONZ site in June 2023 (https://data.neonscience.org/data-products/DP1.00024.001/RELEASE-2023)
    par_data = nu.load_by_product(
        dpid="DP1.00024.001",
        site="KONZ",
        startdate="2023-06",
        enddate="2023-06",
        token=None,
        check_size=False,
        progress=False,
    )

    # Select PAR Data per minute
    par_data = par_data["PARPAR_1min"]
    # Ensure startDateTime is datetime type
    par_data["startDateTime"] = pd.to_datetime(par_data["startDateTime"], utc=True)
    # Convert startDateTime to America/Chicago timezone
    par_data["startDateTime"] = par_data["startDateTime"].dt.tz_convert("America/Chicago").dt.tz_localize(None)
    # Set startDateTime as index and drop columns with all NaN values
    par_data = par_data.set_index("startDateTime").dropna(axis=1, how="all")
    # Locate Date at 19.06
    day_data = par_data.loc["2023-06-19"]
    # Use only the same position values
    day_data = day_data[day_data["horizontalPosition"] == "000"]
    day_data = day_data[day_data["verticalPosition"] == "010"]
    # Limit PPFD to values above 40 to focus on the day period
    day_data = day_data[day_data["PARMean"] >= 40]

    fig, ax = plt.subplots()
    # Plot PAR data
    ax.fill_between(day_data.index, day_data["PARMean"], alpha=0.3, color="black", lw=0)
    # Format axis'
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlabel("Time [hh:mm]")
    ax.set_ylabel(r"PPFD [$\mathrm{\mu mol \, m^{-2} \, s^{-1}}$]")

    # Create day simulation that folows the PAR data
    s = Simulator(model)
    day_prtc = make_protocol(
        [(60, {pfd: row["PARMean"]}) for index, row in day_data.iterrows()]
    )
    res_prior = None
    time_points = 0
    
    while (res_prior is None or isinstance(res_prior, Exception)) and time_points < 1e4:
        s.clear_results()
        time_points += 100
        clear_output(wait=True)
        print(f"Trying simulation with {time_points} time points per step.")
        s.simulate_protocol(day_prtc, time_points_per_step=time_points)
        if isinstance(s.get_result().value, Exception):
            res_prior = None
        else:
            res_prior = "done"
    
    res = s.get_result().unwrap_or_err()
    variables = res.get_variables()
    variables.index = pd.to_datetime(variables.index, unit="s", origin=day_data.index[0])
    fluxes = res.get_fluxes()
    fluxes.index = pd.to_datetime(fluxes.index, unit="s", origin=day_data.index[0])

    res_dict = {}

    for name, pointer in zip(
        [vc, atp, nadph, fluorescence], ["Vc", "ATP", "NADPH", "Fluorescence"]
    ):
        if name is None:
            data = None
            unit = None
        elif name in variables.columns:
            data = variables[name]
            # unit = model._variables[name].unit Reimplement when sperated variables
        elif name in fluxes.columns:
            data = fluxes[name]
            unit = model._reactions[name].unit
        else:
            data = model._parameters[name].value
            unit = model._parameters[name].unit

        res_dict[pointer] = {"data": data, "unit": unit}

    # Colors of sim results
    vc_color = "#fa9442"
    atp_nadph_color = "#008aa1"
    fluo_color = "#DB5764"
    color_list = [vc_color, atp_nadph_color, fluo_color]

    axes_pos = 0.15
    yax_list = []

    #Set lim
    ax.set_xlim(day_data.index[0], day_data.index[-1])
    ax.set_ylim(bottom=0)
    # Create twin axis for each variable to plot
    for ax_idx, color in enumerate(color_list):
        ax_new = ax.twinx()
        ax_new.spines["right"].set_color(color)
        ax_new.spines["right"].set_position(("axes", 1 + axes_pos * ax_idx))
        ax_new.tick_params(axis="y", colors=color)
        yax_list.append(ax_new)

    # Plot variables if they are in the model else write n.a. as ylabel
    vc_ylabel = "Rubisco Carboxylase Activity"
    atp_nadph_ylabel = "ATP/NADPH Ratio"
    fluo_ylabel = "Fluorescence"
    if vc is not None:
        yax_list[0].plot(res_dict["Vc"]["data"], color=vc_color)
        yax_list[0].set_ylabel(
            rf"Rubisco Carboxylase Activity [${custom_latex(res_dict['Vc']['unit'])}$]",
            color=vc_color,
        )
    else:
        vc_ylabel = strike(vc_ylabel)
    yax_list[0].set_ylabel(vc_ylabel, color=vc_color)

    if atp is not None and nadph is not None:
        yax_list[1].plot(
            res_dict["ATP"]["data"] / res_dict["NADPH"]["data"], color=atp_nadph_color
        )
    else:
        atp_nadph_ylabel = strike(atp_nadph_ylabel)
    yax_list[1].set_ylabel(atp_nadph_ylabel, color=atp_nadph_color)

    if fluorescence is not None:
        yax_list[2].plot(res_dict["Fluorescence"]["data"], color=fluo_color)
    else:
        fluo_ylabel = strike(fluo_ylabel)
    yax_list[2].set_ylabel(fluo_ylabel, color=fluo_color)

    return fig, ax