import copy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Model as LmfitModel
from lmfit import Parameters, minimize
from mxlpy import Model, Simulator, make_protocol, plot

from .utils import calc_pam_vals2, create_pamprotocol_from_data, pam_sim

    
def plot_pamfit(
    model: Model,
    new_params: dict | None,
    pfd_str: str,
    fluorescence_str: str | None,
    npq_str: str | None,
    fit_protocol: list[tuple[float, dict]],
    fluo_data: pd.DataFrame,
    sp_lenth: float,
) -> tuple[plt.Figure, plt.Axes]:
    
    F = None
    Fm = None
    NPQ = None
    F_old = None
    Fm_old = None
    NPQ_old = None
    
    if new_params is not None:
        fitted_model = copy.deepcopy(model)
        fitted_model.update_parameters(new_params)
        
        res = pam_sim(
            fit_protocol=make_protocol(fit_protocol),
            model=fitted_model,
            pfd_str=pfd_str,
        )
        # OG model version
        res_old = pam_sim(
            fit_protocol=make_protocol(fit_protocol),
            model=model,
            pfd_str=pfd_str,
        )

        if fluorescence_str is not None:
            F, Fm, NPQ = calc_pam_vals2(res[fluorescence_str], protocol=make_protocol(fit_protocol), pfd_str=pfd_str, do_relative=True)
            F_old, Fm_old, NPQ_old = calc_pam_vals2(res_old[fluorescence_str], protocol=make_protocol(fit_protocol), pfd_str=pfd_str, do_relative=True)
        
        if npq_str is not None:
            NPQ = res[npq_str]
            NPQ_old = res_old[npq_str]
        
    
    data_color = "#84569F"
    fitted_color = "#D64550"
    og_color = "#5B6C5D"
    
    fig, axs = plt.subplot_mosaic([["Fluo", "NPQ"], ["Diff", "Diff"]], figsize=(10, 5))
    
    cleaned_prtc = make_protocol(fit_protocol)
    cleaned_prtc = cleaned_prtc[cleaned_prtc[pfd_str] != 5000]
    
    t_before = 0
    for t, vals in cleaned_prtc.iterrows():
        if vals[pfd_str] == 40:
            color = "#b3b3b3ff"
        elif vals[pfd_str] == 90:
            color = "#d9d9d9ff"
        else:
            color = "#f2f2f2ff"
        for ax in axs:
            axs[ax].set_xlim(0, fluo_data.index[-1])
            axs[ax].axvspan(t_before, t.total_seconds(), facecolor=color, edgecolor="none")
        t_before = t.total_seconds()
    
    #Fitted Data
    axs["Fluo"].plot(F, label="Fitted Fluo", color=fitted_color) if F is not None else None
    axs["Fluo"].plot(Fm, label="Fitted Fm", lw=0, marker="^", color=fitted_color) if Fm is not None else None
    if NPQ is not None:
        if len(NPQ) > len(fluo_data):
            NPQ_scatter = NPQ.loc[fluo_data.index]
        else:
            NPQ_scatter = NPQ
        axs["NPQ"].scatter(NPQ_scatter.index, NPQ_scatter, color=fitted_color, marker="+") if NPQ is not None else None
    axs["NPQ"].plot(NPQ, label="Fitted NPQ", lw=1, color=fitted_color) if NPQ is not None else None
    
    # Experimental Data
    axs["Fluo"].plot(fluo_data.index, fluo_data["F1"], label="Measured Fluo", lw=0, marker="o", color=data_color)
    axs["Fluo"].plot(fluo_data.index + sp_lenth, fluo_data["Fm'1"], label="Measured Fm'", lw=0, marker="x", color=data_color)
    axs["NPQ"].plot(fluo_data["NPQ3"], label="Measured NPQ", lw=1, marker="s", color=data_color)
    data_rel = fluo_data["NPQ3"].copy()
    data_rel = data_rel.replace(0, 1)
    data_rel = data_rel / data_rel
    axs["Diff"].plot(data_rel, label="Measured Baseline", lw=1, ls="dashed", color=data_color)
    
    if NPQ is not None and NPQ_old is not None:
        if len(NPQ) > len(fluo_data):
            NPQ_scatter = NPQ.loc[fluo_data.index]
            NPQ_old_scatter = NPQ_old.loc[fluo_data.index]
        else:
            NPQ_scatter = NPQ
            NPQ_old_scatter = NPQ_old
        axs["Diff"].plot(fluo_data.index, NPQ_scatter / fluo_data["NPQ3"].values, label="Fitted / Measured", lw=1, marker="o", color=fitted_color)
        axs["Diff"].plot(fluo_data.index, NPQ_old_scatter / fluo_data["NPQ3"].values, label="OG Model / Measured", lw=1, marker="x", color=og_color)
    
    axs["Diff"].set_ylim(-0.1, 2.1)
    axs["Diff"].set_ylabel("Model / Measured NPQ [a.u.]")
    axs["Diff"].set_xlabel("Time [s]")
    
    axs["Fluo"].set_xlabel("Time [s]")
    axs["Fluo"].set_ylabel("Fluorescence [a.u.]")
    
    
    axs["NPQ"].set_xlabel("Time [s]")
    axs["NPQ"].yaxis.tick_right()
    axs["NPQ"].yaxis.set_label_position("right")
    axs["NPQ"].set_ylabel("NPQ [a.u.]")
    
    axs["Fluo"].legend(loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1))
    axs["NPQ"].legend(loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1))
    axs["Diff"].legend(loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1))
    
    param_str = "Relative Differences:\n"
    
    if new_params is not None:
        for param, val in new_params.items():
            rel_val = ((val - model.get_raw_parameters()[param].value) / model.get_raw_parameters()[param].value) * 100
            rel_val = f"+{rel_val:.2f}%" if rel_val >= 0 else f"{rel_val:.2f}%"
            param_str += f"{param} = {rel_val}\n"
        
    plt.text(0.9, 0.5, param_str, transform=fig.transFigure, ha="left", va="center")
    
    plt.tight_layout()
    
    return fig, axs

def pamfit_lmmodel(
    npq_x: pd.Series,
    npq_y_mean: pd.Series,
    model: Model,
    fit_protocol_dict: dict,
    pfd_str: str,
    fluorescence_str: str | None,
    npq_str: str | None,
    relative: bool,
    sat_pulse: float,
    **params
) -> np.ndarray:
    fit_protocol = pd.DataFrame(fit_protocol_dict)
    fit_protocol.index.name = "Timedelta"
    model = copy.deepcopy(model)
    model.update_parameters(params)
    
    res = pam_sim(
        fit_protocol=fit_protocol,
        model=model,
        pfd_str=pfd_str,
    )
    
    if res is None:
        return np.ones(len(npq_x)) * 1e6
    
    if npq_str is not None:
        NPQ = res[npq_str].loc[npq_x]
    else:
        F, Fm, NPQ = calc_pam_vals2(
                fluo_result=res[fluorescence_str],
                protocol=fit_protocol,
                pfd_str=pfd_str,
                sat_pulse=sat_pulse,
                do_relative=relative
            )
        
    return NPQ.values - npq_y_mean

def create_pamfit(
    model: Model,
    pfd_str: str,
    fluorescence_str: str | None,
    npq_str: str | None,
    pam_params_to_fit: list[str],
    relative: bool = True,
    standard_scale: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    fluo_data = pd.read_csv(Path(__file__).parent / "Data/fluo_col0_1.csv", index_col=0) # Taken from https://doi.org/10.1111/nph.18534
    # Data taken with Maxi Imaging-PAM (Walz, Germany) using Col-0 Arabidopsis thaliana plants.
    # SP standard length = 720ms Maximal setting is standard (level 10) = 5000 µmol m-2 s-1 on IMAG-MAX/L TODO: Check if it is correct
    fluo_data["F1"] = fluo_data["F1"] / fluo_data["Fm'1"].iloc[0]
    fluo_data["Fm'1"] = fluo_data["Fm'1"] / fluo_data["Fm'1"].iloc[0]
    fluo_data["NPQ3"] = (fluo_data["Fm'1"].iloc[0] - fluo_data["Fm'1"]) / fluo_data["Fm'1"]
    
    sp_lenth = 720 / 1000  # seconds
    sp_intensity = 5000  # µmol m-2 s-1
    
    #Convert index to time in seconds
    fluo_data.index = pd.to_timedelta(fluo_data.index)
    fluo_data.index = fluo_data.index - fluo_data.index[0]
    fluo_data.index = fluo_data.index.total_seconds()
    
    # Complete standard scaling if required
    data_mean = fluo_data["NPQ3"].mean() if standard_scale else 0
    data_std = fluo_data["NPQ3"].std() if standard_scale else 1
    fluo_data["NPQ_central"] = (fluo_data["NPQ3"] - data_mean) if standard_scale else fluo_data["NPQ3"]
    
    # Create Pam Protocol
    fit_protocol = create_pamprotocol_from_data(
        data=fluo_data,
        par_column="PAR",
        pfd_str=pfd_str,
        time_sp=sp_lenth, #720ms SP to seconds
        sp_pluse=sp_intensity # 5000 µmol m-2 s-1 on IMAG-MAX/L
    )


    if fluorescence_str is not None or npq_str is not None:
        fit_model = LmfitModel(
            func=pamfit_lmmodel,
            independent_vars=["npq_x", "npq_y_mean", "model", "fit_protocol_dict", "pfd_str", "fluorescence_str", "npq_str", "relative", "sat_pulse"],
        )
    
        initial_params = Parameters()
        for param in pam_params_to_fit:
            val = model.get_raw_parameters()[param].value
            if val < 0:
                initial_params.add(param, value=val, vary=True, max=0)
            else:
                initial_params.add(param, value=val, vary=True, min=0)
    
        result = fit_model.fit(
            data=fluo_data["NPQ_central"].values,
            params=initial_params,
            weights=1 / data_std,
            npq_x=fluo_data["NPQ_central"].index,
            npq_y_mean=data_mean,
            model=model,
            fit_protocol_dict=make_protocol(fit_protocol).to_dict(),
            pfd_str=pfd_str,
            fluorescence_str=fluorescence_str,
            npq_str=npq_str,
            relative=relative,
            sat_pulse=sp_intensity,
        )
        best_params = result.best_values
        print(best_params)
    else:
        best_params = None
    
    fig, axs = plot_pamfit(
        model=model,
        new_params=best_params,
        pfd_str=pfd_str,
        fluorescence_str=fluorescence_str,
        npq_str=npq_str,
        fit_protocol=fit_protocol,
        fluo_data=fluo_data,
        sp_lenth=sp_lenth,
    )

    return fig, axs