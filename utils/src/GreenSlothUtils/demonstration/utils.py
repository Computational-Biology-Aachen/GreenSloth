import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import Normalize
from mdutils.mdutils import MdUtils
from mxlpy import Model, Simulator, mca, plot
from scipy.signal import find_peaks, peak_prominences
from sympy import Integer, Pow
from sympy.physics.units import Unit
from sympy.printing.latex import LatexPrinter
from IPython.display import clear_output


# Custom Latex Printer to handle units with negative integer exponents properly (Still need to be improved for more complex cases)
class RefinedUnitLatexPrinter(LatexPrinter):
    def _print_Pow(self, expr):
        # Custom handling for Unit objects with negative integer exponents
        if (
            isinstance(expr.base, Unit)
            and isinstance(expr.exp, Integer)
            and expr.exp < 0
        ):
            base_str = self._print(expr.base)
            exp_str = self._print(expr.exp)
            return f"{base_str}^{{{exp_str}}}"
        # Fallback to default behavior for other Pow expressions
        return super()._print_Pow(expr)

    def _print_Mul(self, expr):
        # Separate terms with positive and negative exponents
        positive_terms = []
        negative_terms = []

        for arg in expr.args:
            # Check if the argument is a Pow of a Unit with a negative integer exponent
            if (
                isinstance(arg, Pow)
                and isinstance(arg.base, Unit)
                and isinstance(arg.exp, Integer)
                and arg.exp < 0
            ):
                negative_terms.append(arg)
            else:
                positive_terms.append(arg)

        # Print positive terms and negative terms
        # SymPy's default _print(Unit) for millimeter is \text{mm}, so no extra \cdot is needed within it.
        printed_positive_terms = [self._print(t) for t in positive_terms]
        printed_negative_terms = [self._print(t) for t in negative_terms]

        all_printed_terms = []
        if printed_positive_terms:
            all_printed_terms.extend(printed_positive_terms)
        if printed_negative_terms:
            all_printed_terms.extend(printed_negative_terms)

        if all_printed_terms:
            return r"\ ".join(all_printed_terms)
        else:
            # Fallback to default behavior if no specific handling is needed
            return super()._print_Mul(expr)


def custom_latex(expr):
    return RefinedUnitLatexPrinter().doprint(expr)


def calc_co2_conc(pco2: float, H_cp_co2: float = 3.4e-4):
    """Calculate the CO2 concentration based on CO2 partial pressure and Henry's law constant.

    Args:
        pco2 (float): CO2 partial pressure [µbar].
        H_cp_co2 (float, optional): Henry's law constant for CO2 at 25°C [mM * Pa-1]. Defaults to 3.4e-4 (https://doi.org/10.5194/acp-23-10901-2023).

    Returns:
        float: CO2 concentration in mM.
    """
    # Unit conversions
    H_cp_co2 = H_cp_co2 * 1e5  # [mM * bar-1]
    H_cp_co2 = H_cp_co2 * 1e-6  # [mM * µbar-1]

    return H_cp_co2 * pco2


def mM_to_µmol_per_m2(conc_mM: float, corr_factor: float = 0.0112):
    """Convert mM concentration to µmol m-2.

    Args:
        conc_mM (float): Concentration in mM.
        corr_factor (float, optional): Correction factor. Defaults to 0.0112, which is the factor for the stroma (https://doi.org/10.1007/s11120-006-9109-1).

    Returns:
        float: Concentration in µmol m-2.
    """
    return conc_mM * 1e3 * corr_factor

def calc_pam_vals(
    fluo_result: pd.Series, peak_distance: float = 120
) -> tuple[pd.Series, pd.Series, pd.Series, pd.DataFrame]:
    """Calculate PAM values from fluorescence data.

    Use the fluorescence data from a PAM protocol to calculate Fm, NPQ, Fmin, and the quantum yields Y(NO), Y(NPQ) and Y(II). To find the Fm values, the peaks in the fluorescence data are found using scipy.signal.find_peaks. The distance between the peaks should be the same length as a period used in the PAM protocol, however may need to be adjusted based on the fluorescence data. Best to plot the Fluorescence data and the calculated Fm to check if the peaks are found correctly.

    Calculations and Assumptions:
        Fm(t): Maximum fluorescence at time t
        Fmin(t): Minimum fluorescence left of peak at time t
        NPQ(t) = (Fm(0) - Fm(t)) / Fm(t): Non-photochemical quenching at time t (Fm(0) is the first Fm value and is assumed to be the highest Fm value)
        Y(NO)(t) = Fmin(t) / Fm(0) : Quantum Yield of Non-Regulated Energy Loss at time t
        Y(NPQ)(t) = Fmin(t) / Fm(t) - Fmin(t) / Fm(0) : Quantum Yield of Regulated Heat Dissipation at time t
        Y(II)(t) = (Fm(t) - Fmin(t)) / Fm(t) : Quantum Yield of Photochemical Energy Conversion at time t

    Args:
        fluo_result (pd.Series): Fluorescence data as a pd.Series from mxlpy simulation.
        peak_distance (float, optional): Minimum distance between peaks, which should be the same length as a period used in the PAM protocol. However may need to be adjusted based on the fluorescence data. Defaults to 120 s.

    Returns:
        Fm (pd.Series): Maximum fluorescence values
        NPQ (pd.Series): Non-photochemical quenching values
        Fmin (pd.Series): Minimum fluorescence values
        quant_yields (pd.DataFrame): Quantum yields (Y(NO), Y(NPQ), Y(II))
    """

    # Find the indices of the Flourescence peaks (Fmaxs)
    peaks, _ = find_peaks(fluo_result, distance=peak_distance, height=0)

    # Fm series
    Fm = fluo_result.iloc[peaks]
    Fm.name = "Flourescence Peaks (Fm)"

    # Calculate NPQ
    NPQ = (Fm.iloc[0] - Fm) / Fm if len(Fm) > 0 else pd.Series(dtype=float)
    NPQ.name = "Non-Photochemical Quenching (NPQ)"

    # Find the minima around the peaks
    prominences, prominences_left, prominences_right = peak_prominences(
        (fluo_result), peaks, wlen=peak_distance
    )

    # Fmin is always the minima before the peak
    Fmin = fluo_result.iloc[prominences_left]
    Fmin.name = "Minimum Fluorescence before peaks (Fmin)"

    # Quantum Yield of Non-Regulated Energy Loss (Y(NO))
    Y_NO = Fmin / Fm.iloc[0] if len(Fm) > 0 else pd.Series(dtype=float)
    Y_NO.name = "Y(NO)"

    # Quantum Yield of Regulated Heat Dissipation (Y(NPQ))
    Y_NPQ = (
        Fmin / Fm.values - Fmin / Fm.iloc[0] if len(Fm) > 0 else pd.Series(dtype=float)
    )
    Y_NPQ.name = "Y(NPQ)"

    # Quantum Yield of Photochemical Energy Conversion (Y(II))
    Y_II = (Fm.values - Fmin) / Fm.values if len(Fm) > 0 else pd.Series(dtype=float)
    Y_II.name = "Y(II)"

    # pd.DataFrame of the three quantum yields
    quant_yields = pd.concat([Y_NO, Y_NPQ, Y_II], axis=1)
    quant_yields["Total"] = quant_yields.sum(axis=1)

    return Fm, NPQ, Fmin, quant_yields

def calc_pam_vals2(
    fluo_result: pd.Series, protocol: pd.DataFrame, pfd_str: str, sat_pulse: float = 2000, do_relative: bool = False
) -> tuple[pd.Series, pd.Series]:
    """Calculate PAM values from fluorescence data.

    Use the fluorescence data from a PAM protocol to calculate Fm, NPQ. To find the Fm values, the protocol used for simulation is seperated into ranges between each saturating pulse. Then the maximum fluorescence value within each range is taken as Fm. Thes are then used to calculate NPQ.

    Args:
        fluo_result (pd.Series): Fluorescence data as a pd.Series from mxlpy simulation.
        protocol (pd.DataFrame): PAM protocol used for simulation. Created using make_protocol from mxlpy.
        pfd_str (str): The name of the PPFD parameter in the protocol.
        sat_pulse (float, optional): The threshold for saturating pulse in the protocol. Defaults to 2000.

    Returns:
        tuple[pd.Series, pd.Series]: Fm and NPQ as pd.Series
    """    
    
    F = fluo_result.copy()
    F.name = "Fluorescence"
    
    peaks = protocol[protocol[pfd_str] >= sat_pulse].copy()
    peaks.index = peaks.index.total_seconds()
    peaks = peaks.reset_index()
    
    Fm = {
        "start": [],
        "end": [],
        "time": [],
        "value": []
    }

    for idx, (time, _) in peaks.iterrows():
        if idx == 0:
            start_time = 0
        else:
            start_time = time - (time - peaks["Timedelta"].iloc[idx - 1]) / 2
            
        if idx == len(peaks) - 1:
            end_time = fluo_result.index[-1]
        else:
            end_time = time + (peaks["Timedelta"].iloc[idx + 1] - time) / 2
            
        Fm["start"].append(start_time)
        Fm["end"].append(end_time)
        Fm_slice = fluo_result.loc[start_time:end_time]
        Fm["time"].append(Fm_slice.idxmax())
        Fm["value"].append(Fm_slice.max())
        
    Fm = pd.DataFrame(Fm).set_index("time")
    Fm = Fm["value"]
    Fm.name = "Flourescence Peaks (Fm)"
    
    if do_relative:
        F = F / Fm.iloc[0]
        Fm = Fm / Fm.iloc[0]
    
    # Calculate NPQ
    NPQ = (Fm.iloc[0] - Fm) / Fm if len(Fm) > 0 else pd.Series(dtype=float)
    NPQ.name = "Non-Photochemical Quenching (NPQ)"
    
    return F, Fm, NPQ


def create_pamprotocol_from_data(
    data: pd.DataFrame,
    par_column: str,
    pfd_str: str,
    time_sp: float,
    sp_pluse: float
):
    time_simed = 0
    fit_protocol = []
    dark_light = 40
    
    for time in data.index:
        if time != 0:
            if data.loc[time, par_column] == 0:
                pfd_val = dark_light
            else:
                pfd_val = data.loc[time, par_column]
            fit_protocol.append((time - time_simed - time_sp, {pfd_str: pfd_val}))
        fit_protocol.append((time_sp, {pfd_str: sp_pluse}))
        time_simed = time
        
    return fit_protocol

def pam_sim(
    fit_protocol: list[tuple[float, dict]],
    model: Model,
    pfd_str: str,
    dark_adaptation_time: float = 60*30,
    dark_pfd: float = 40,
) -> pd.DataFrame:
    s = Simulator(model=model)
    
    s.update_parameter(pfd_str, dark_pfd)
    res_prior = None
    time_points = 0
    
    while (res_prior is None or isinstance(res_prior, Exception)) and time_points < 1e4:
        s.clear_results()
        time_points += 100
        clear_output(wait=True)
        print(f"Trying dark-simulation with {time_points} time points per step.")
        s.simulate(dark_adaptation_time, steps=time_points)
        if isinstance(s.get_result().value, Exception):
            res_prior = None
        else:
            res_prior = "done"
        
    if res_prior is None:
        print("No result from dark simulation")
        return None
    
    
    dark_y0 = s.get_result().unwrap_or_err().get_new_y0()
    
    s.clear_results()
    s.update_variables(dark_y0)
        
    res = None
    time_points = 0
    
    while (res is None or isinstance(res, Exception)) and time_points < 1e4:
        s.clear_results()
        time_points += 100
        clear_output(wait=True)
        print(f"Trying pam-simulation with {time_points} time points per step.")
        s.simulate_protocol(fit_protocol, time_points_per_step=time_points)
        res = None if isinstance(s.get_result().value, Exception) else "done"

    res = s.get_result().unwrap_or_err().get_combined()
    
    return res

def param_recursion(
    model: Model,
    search_str: str,
    dict_out: dict,
    order: int = 0,
    max_order: int = 5,
):
    """Recursion fucntion to find parameters influencing a given variable, reaction, derived variable or readout in a MxLpy model.

    Args:
        model (Model): mxlpy model to recursivle search for parameters.
        search_str (str): Entity to search for parameters influencing it. Needs to be in the model.
        dict_out (dict): Dictionary to store found parameters in. Should already exist when calling the function.
        order (int, optional): Number of recursion order. Defaults to 0.
        max_order (int, optional): Maximum recursion order. Defaults to 5.
    """    
    if order > max_order:
        return
    dict_key = f"Order {order}"
    if dict_out.get(dict_key) is None:
        dict_out[f"Order {order}"] = []
        
    type_of_id = model.ids[search_str]

    if type_of_id == "readout":
        to_fit = model._readouts[search_str]
    elif type_of_id == "variable":
        stoics = model.get_stoichiometries_of_variable(search_str)
        for reac in stoics.keys():
            param_recursion(model, reac, order=order+1, dict_out=dict_out, max_order=max_order)
        return
    elif type_of_id == "derived":
        to_fit = model._derived[search_str]
    elif type_of_id == "reaction":
        to_fit = model._reactions[search_str]

    for arg in to_fit.args:
        if arg not in model.ids:
            continue
        elif model.ids[arg] == "parameter":
            dict_out[f"Order {order}"].append(arg)
        else:
            param_recursion(model, arg, order=order+1, dict_out=dict_out, max_order=max_order)

def find_params_to_fit_byorder(
    to_fit_str: str,
    model: Model,
    max_order: int = 5,
) -> None:
    """Recursviely looks through provided model to find parameters influencing a given variable, reaction, derived variable or readout. It will print the parameters found at each order up to the given maximum order. The smaller the order the nearer the parameter is to the fitted entity.

    Args:
        to_fit_str (str): Name of variable, reaction, derived variable or readout to fit. It needs to be in the model.
        model (Model): The mxlpy model that should be anlyzed fro fitting.
        max_order (int, optional): The maximum order of parameters to find. Defaults to 5.
    """    
    dict_out = {}
    param_recursion(model, to_fit_str, order=1, dict_out=dict_out, max_order=max_order)
    max_order_length = max([len(v) for v in dict_out.values()])
    
    vars_rcoeffs, flux_rcoeffs = mca.response_coefficients(
        model=model,
        to_scan=None,
        normalized=True,
    )
    
    for rcoeffs in [vars_rcoeffs, flux_rcoeffs]:
        if to_fit_str not in rcoeffs.index:
            continue
        
        correct_coeffs = rcoeffs.loc[to_fit_str]
    
    fig, axs = plt.subplot_mosaic(
        mosaic=[[f"Order {i+1}", "cbar"] for i in range(max_order)],
        width_ratios=[1, 0.2],
        figsize=(6, max_order * 4)
    )
    
    for i in range(max_order):
        order_str = f"Order {i+1}"
        lst_params = dict_out.get(order_str, [])
        plot.heatmap(
            ax=axs[order_str],
            df=correct_coeffs[lst_params].to_frame(),
            invert_yaxis=False,
            annotate=True,
            norm=Normalize(vmin=-1, vmax=1),
            cax=axs["cbar"] if i == 0 else None,
            colorbar=False if i != 0 else True,
        )
        axs[order_str].set_xticks([])
        axs[order_str].set_title(order_str)
        
    plt.tight_layout()
    
    plt.show()
    
def find_params_to_fit_byelasticities(
    to_fit_str: str,
    model: Model,
    max_num: int = 30,
    omit_strs: None | list[str] = None,
):   
    # dict_out = {}
    # param_recursion(model, to_fit_str, order=1, dict_out=dict_out, max_order=max_order)
    # max_order_length = max([len(v) for v in dict_out.values()])
    
    vars_rcoeffs, flux_rcoeffs = mca.response_coefficients(
        model=model,
        to_scan=None,
        normalized=True,
    )
    
    for rcoeffs in [vars_rcoeffs, flux_rcoeffs]:
        if to_fit_str not in rcoeffs.index:
            continue
        
        correct_coeffs = rcoeffs.loc[to_fit_str].to_frame()
        
    correct_coeffs["Abs"] = correct_coeffs[to_fit_str].abs()
    correct_coeffs = correct_coeffs.sort_values(by="Abs", ascending=False)
        
    sorted_coeffs = correct_coeffs[to_fit_str]
    if omit_strs is not None:
        sorted_coeffs = sorted_coeffs[[i for i in sorted_coeffs.index if i not in omit_strs]]
    sorted_coeffs = sorted_coeffs.iloc[:max_num]
    
    fig, ax, _ = plot.heatmap(
        df=sorted_coeffs.to_frame(),
        invert_yaxis=True,
        annotate=True,
    )
    
    ax.set_xticks([])
    # axs[order_str].set_title(order_str)
        
    plt.tight_layout()
    
    plt.show()

def save_matplotlib_figure(fig: plt.Figure, file_prepend: str, figcat: str) -> None:
    """Save matplotlib figure of model validation fig in 'Figures' directory as a svg.

    Args:
        fig (plt.Figure): Matplotlib figure to save
        file_prepend (str): str of file prepend for each fig. Recommended to be the model name.
        figcat (str): str of figure category. Recommended to be a descriptive name for the figure.
    """

    try:
        notebook_dir = Path(os.getcwd())
    except NameError:
        notebook_dir = Path(".")

    # 2. Construct the path for the 'Figure' directory
    figure_dir = notebook_dir / "Figures"

    # 3. Create the directory if it doesn't exist
    # `mkdir(exist_ok=True)` prevents an error if the directory already exists.
    figure_dir.mkdir(parents=True, exist_ok=True)

    # 4. Construct the full save path
    save_path = figure_dir / f"{file_prepend}_{figcat}.svg"

    # 5. Save the figure
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)  # Close the figure to free up memory


def create_save_figs(
    model: Model,
    pfd: str,
    file_prepend: str,
    co2: str | None,
    vc: str | None,
    pco2: str | None,
    H_cp_co2: str | None,
    gammastar: str | None,
    r_light: str | None,
    A: str | None,
    flourescence: str | None,
    atp: str | None,
    nadph: str | None,
    rubp: str | None,
    coeff_psii: str | None,
    coeff_psi: str | None,
    coeff_rubisco: str | None,
) -> None:
    """Create and save all model validation figs.

    Use a MxLpy model to create and save all model validation figs. The figs include:
    - FvCB comparison fig
    - PAM simulation fig
    - Day simulation fig
    - MCA of photosynthesis control coefficients fig

    Args:
        model (Model): MxLpy model to create figs from.
        pfd (str): Name for PPFD parameter in model.
        file_prepend (str): str of file prepend for each fig. Recommended to be the model name.
        co2 (str | None): Name for CO2 variable in model.
        vc (str | None): Name for rubisco carboxylation in model.
        pco2 (str | None): Name for CO2 partial pressure in model.
        H_cp_co2 (str | None): Name for Henry's law constant for CO2 in model.
        gammastar (str | None): Name for CO2 compensation point in model.
        r_light (str | None): Name for Rate of non-photorespiratory CO2 release in the light in model.
        A (str | None): Name for net carbon assimilation rate in model.
        flourescence (str | None): Name for fluorescence in model.
        atp (str | None): Name for ATP in model.
        nadph (str | None): Name for NADPH in model.
        rubp (str | None): Name for RuBP in model.
        coeff_psii (str | None): Name for PSII coefficient in model.
        coeff_psi (str | None): Name for PSI coefficient in model.
        coeff_rubisco (str | None): Name for Rubisco coefficient in model.
    """

    # FvCB comparison
    plot, ax = create_fvcb_fig(
        model=model,
        pfd=pfd,
        co2=co2,
        vc=vc,
        pco2=pco2,
        H_cp_co2=H_cp_co2,
        gammastar=gammastar,
        r_light=r_light,
        A=A,
    )
    save_matplotlib_figure(plot, file_prepend, "fvcb_compare")

    # PAM simulation
    pam_plot, ax = create_pam_fig(model=model, pfd=pfd, flourescence=flourescence)
    save_matplotlib_figure(pam_plot, file_prepend, "pam")

    # Day simulation
    plot, ax = create_day_simulation_fig(
        model=model, pfd=pfd, vc=vc, atp=atp, nadph=nadph, flourescence=flourescence
    )
    save_matplotlib_figure(plot, file_prepend, "day_simulation")

    # MCA of photosynthesis control coefficients
    plot, ax = create_mca_fig(
        model=model,
        coeff_psii=coeff_psii,
        coeff_psi=coeff_psi,
        coeff_rubisco=coeff_rubisco,
        rubp=rubp,
        co2=co2,
    )

    return


def create_report_summary(
    model: Model,
    pfd: str,
    file_prepend: str,
    co2: str | None = None,
    vc: str | None = None,
    pco2: str | None = None,
    H_cp_co2: str | None = None,
    gammastar: str | None = None,
    r_light: str | None = None,
    A: str | None = None,
    flourescence: str | None = None,
    atp: str | None = None,
    nadph: str | None = None,
    rubp: str | None = None,
    coeff_psii: str | None = None,
    coeff_psi: str | None = None,
    coeff_rubisco: str | None = None,
):
    """Create markdown file of Model validation.

    Use a MxLpy model to create and save all model validation figs and create a markdown summary file. The figs include:
    - FvCB comparison fig
    - PAM simulation fig

    Args:
        model (Model): MxLpy model to create figs from.
        pfd (str): Name for PPFD parameter in model.
        file_prepend (str): str of file prepend for each fig. Recommended to be the model name.
        co2 (str | None, optional): Name for CO2 variable in model. Defaults to None.
        vc (str | None, optional): Name for rubisco carboxylation in model. Defaults to None.
        pco2 (str | None, optional): Name for CO2 partial pressure in model. Defaults to None.
        H_cp_co2 (str | None, optional): Name for Henry's law constant for CO2 in model.. Defaults to None.
        gammastar (str | None, optional): Name for CO2 compensation point in model. Defaults to None.
        r_light (str | None, optional): Name for Rate of non-photorespiratory CO2 release in the light in model. Defaults to None.
        A (str | None, optional): Name for net carbon assimilation rate in model. Defaults to None.
        flourescence (str | None, optional): Name for fluorescence in model. Defaults to None.
        atp (str | None, optional): Name for ATP in model. Defaults to None.
        nadph (str | None, optional): Name for NADPH in model. Defaults to None.
        rubp (str | None, optional): Name for RuBP in model. Defaults to None.
        coeff_psii (str | None, optional): Name for PSII coefficient in model. Defaults to None.
        coeff_psi (str | None, optional): Name for PSI coefficient in model. Defaults to None.
        coeff_rubisco (str | None, optional): Name for Rubisco coefficient in model. Defaults to None.
    """

    # TODO: Add PAM assumptions
    # TODO: Add Day simulation
    # TODO: Add MCA

    # Create all the validation figs and save them
    create_save_figs(
        model=model,
        pfd=pfd,
        file_prepend=file_prepend,
        co2=co2,
        vc=vc,
        pco2=pco2,
        H_cp_co2=H_cp_co2,
        gammastar=gammastar,
        r_light=r_light,
        A=A,
        flourescence=flourescence,
        atp=atp,
        nadph=nadph,
        rubp=rubp,
        coeff_psii=coeff_psii,
        coeff_psi=coeff_psi,
        coeff_rubisco=coeff_rubisco,
    )

    # Create markdown summary file
    mdFile = MdUtils(
        file_name=f"{Path(__file__).parent / f'{file_prepend}_report_summary.md'}",
        title=f"{file_prepend} Report Summary",
    )
    mdFile.new_header(level=1, title="Simulations")

    # Carbon Assimilation via FvCB
    mdFile.new_header(level=2, title="Carbon Assimilation via FvCB")

    # Create table of parameters that are in the model
    table = ["Parameter", "Exists?"]

    for text, param in zip(
        [
            r"$\mathrm{CO}_2$",
            r"$v_\mathrm{c}$",
            r"$\mathrm{C_i}$",
            r"$H_\mathrm{s}^{cp}$",
            r"$\Gamma ^*$",
            r"$R_\mathrm{light}$",
            r"A",
        ],
        [co2, vc, pco2, H_cp_co2, gammastar, r_light, A],
    ):
        table.append(text)
        if param is None:
            if text in ["CO2", "Vc"]:
                table.append("&#10060;")
            else:
                table.append("&#10060;")
        else:
            table.append("&check;")

    # Explanation of FvCB comparison
    mdFile.new_paragraph(
        r"Comparison of modelled carbon assimilation ($A$) and carboxylation rate ($v_\mathrm{c}$) against the Farquhar, von Caemmerer and Berry (FvCB) model. The FvCB model is calculated using the min-W approach as described by Lochoki and McGrath (2025) [[1]](https://doi.org/10.1101/2025.03.11.642611). To be able to simulate carbon assimilation, there are two mandatory parameters that need to be present in the model: CO2 concentration and Vc. If one of these parameters is missing, the FvCB model will still be shown, but no comparison with the model will be possible. Other parameters that are required to calculate the FvCB model will be added as parameters with default values if they are not present in the model. The table below summarizes which parameters were found in the model. The carbon assimilation shown does not represent actual values but rather a theoretical curve to compare the kinetic model to the popular FvCB model."
    )

    mdFile.new_header(level=3, title="Assumptions")

    # List of assumptions made for FvCB comparison
    mdFile.new_list(
        [
            r"Infinite mesophyll conductance, therefore intercellular CO<sub>2</sub> partial pressure equals chloroplast partial pressure ($\mathrm{C_i} = \mathrm{C_c}$)",
            r"If no CO<sub>2</sub> concentration nor rate of rubisco carboxylation ($v_\mathrm{c}$) is present in the model, no comparison will be shown",
            r"If no $\mathrm{C_i}$ is present in the model, it will be added as a parameter assuming an initial value of CO<sub>2</sub> concentration divided by Henry's law constant for CO<sub>2</sub> ($H_\mathrm{s}^{cp}$)",
            r"If no $H_\mathrm{s}^{cp}$ is present in the model, it will be added as a parameter with a value of $3.4 \times 10^{-4}\ \mathrm{mmol\ Pa^ {-1}}$ [[2]](https://doi.org/10.5194/acp-23-10901-2023)",
            r"If no CO<sub>2</sub> compensation point in the absence of non-photorespiratory CO<sub>2</sub> release ($\Gamma ^*$) is present in the model, it will be added as a parameter with a value of $38.6\ \mathrm{\mu bar}$ [[1]](https://doi.org/10.1101/2025.03.11.642611)",
            r"If no $R_\mathrm{light}$ is present in the model, it will be added as a parameter with a value of $1\ \mathrm{\mu mol\ m^{-2}\ s^{-1}}$ [[1]](https://doi.org/10.1101/2025.03.11.642611)",
            r"If no $A$ is present in the model, it will be added as a derived variable following the FvCB equation [[1]](https://doi.org/10.1101/2025.03.11.642611): $v_\mathrm{c} \cdot \left(1 - \frac{\Gamma ^*}{C_i}\right) - R_\mathrm{light}$",
            r"To be able to compare with original FvCB curves, the model needs to have $v_\mathrm{c}$ following the same units as the FvCB model ($\mathrm{\mu mol\ m^{-2}\ s^{-1}}$). The `mM_to_µmol_per_m2` can be used to convert from mM to $\mathrm{\mu mol\ m^{-2}}$ assuming a volume factor of $0.0112\ \mathrm{L\ m^{-2}}$ in the stroma [[3]](https://doi.org/10.1007/s11120-006-9109-1). If the given units are in mM, the conversion will be done automatically, by adding a derived parameter with the converted values.",
        ]
    )

    # Add FvCB comparison fig
    mdFile.new_line(
        f"{mdFile.new_inline_image(text='Assimilation', path=str(f'./Figures/{file_prepend}_fvcb_compare.svg'))}"
    )

    mdFile.new_table(columns=2, rows=8, text=table)

    # PAM Fluorescence
    mdFile.new_header(level=2, title="PAM Fluorescence")

    # Explanation of PAM simulation
    mdFile.new_paragraph(
        r"Simulation of a PAM flourescence protocol. The simulation is first run for 30 minutes in a dark adapted state (PPFD = 40) and then the PAM protocol starts. Each period consists of 2 minutes of light and then a saturating pulse of 0.8 seconds. The first two periods are in low light (PPFD = 40), followed by 10 periods in actinic light (PPFD = 1000) and then 10 periods in low light again (PPFD = 40). The left plot shows the normalised flourescence yield (orange) with the identified Fm peaks (crosses) and the calculated NPQ (blue). The right plot shows the quantum yields of non-regulated energy loss (Y(NO), red), regulated heat dissipation (Y(NPQ), orange) and photochemical energy conversion (Y(II), green), but only during the light phase. All results here are arbituary by using the proposed initial conditions of the model and using the Flourescence readout calculated through the model. Therefore, the values do not represent actual values but rather a qualitative behaviour of the model. The table below summarizes which parameters were found in the model."
    )

    # Add Simulation fig
    mdFile.new_line(
        f"{mdFile.new_inline_image(text='PAM Protocol', path=str(f'./Figures/{file_prepend}_pam.svg'))}"
    )

    # Day Simulation
    mdFile.new_header(level=2, title="Day Simulation")

    # Add to markdown file
    mdFile.create_md_file()

    return