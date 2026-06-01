from __future__ import annotations

import ast
import inspect
from collections.abc import Callable
from pathlib import Path
from functools import partial

import mxlbricks
import mxlbricks.names as n
from mxlpy import Model
from mxlpy.types import Derived, InitialAssignment
import pandas as pd

from typing import Literal

EMPTY: Literal[""] = ""


def new_name_creator(new_name: str) -> Callable:
    def new_func(compartment: str = EMPTY, tissue: str = EMPTY) -> str:
        return n.loc(new_name, compartment, tissue)
    return new_func
    

def redefine_names(
    vars_glossary_path: Path,
    rates_glossary_path: Path,
) -> None:
    
    n.pfd = new_name_creator("PPFD")
    
    vars_gloss = pd.read_csv(vars_glossary_path, keep_default_na=False, dtype={"Glossary ID": "Int64"}, usecols=["Python Var", "Glossary ID"])
    
    vars_gloss = vars_gloss.set_index("Glossary ID")
    
    vars_gloss["New Name"] = vars_gloss["Python Var"].apply(new_name_creator)

    n.pq_red = vars_gloss["New Name"].loc[1]
    n.atp = vars_gloss["New Name"].loc[2]
    vars_gloss["New Name"].loc[3] # Not set, as it is based on compartment
    n.psbs_de = vars_gloss["New Name"].loc[4]
    n.vx = vars_gloss["New Name"].loc[5]
    vars_gloss["New Name"].loc[6] # ATPase
    n.pq_ox = vars_gloss["New Name"].loc[7]
    n.adp = vars_gloss["New Name"].loc[8]
    n.psbs_pr = vars_gloss["New Name"].loc[9]
    n.zx = vars_gloss["New Name"].loc[10]
    vars_gloss["New Name"].loc[11] # Atpase_inac
    n.b0 = vars_gloss["New Name"].loc[12]
    n.b1 = vars_gloss["New Name"].loc[13]
    n.b2 = vars_gloss["New Name"].loc[14]
    n.b3 = vars_gloss["New Name"].loc[15]
    vars_gloss["New Name"].loc[16] #pH Lumen
    n.pc_ox = vars_gloss["New Name"].loc[17]
    n.fd_ox = vars_gloss["New Name"].loc[18]
    n.nadph = vars_gloss["New Name"].loc[19]
    n.lhc = vars_gloss["New Name"].loc[20]
    n.pga = vars_gloss["New Name"].loc[21]
    n.bpga = vars_gloss["New Name"].loc[22]
    n.gap = vars_gloss["New Name"].loc[23]
    n.dhap = vars_gloss["New Name"].loc[24]
    n.fbp = vars_gloss["New Name"].loc[25]
    n.f6p = vars_gloss["New Name"].loc[26]
    n.g6p = vars_gloss["New Name"].loc[27]
    n.g1p = vars_gloss["New Name"].loc[28]
    n.sbp = vars_gloss["New Name"].loc[29]
    n.s7p = vars_gloss["New Name"].loc[30]
    n.e4p = vars_gloss["New Name"].loc[31]
    n.x5p = vars_gloss["New Name"].loc[32]
    n.r5p = vars_gloss["New Name"].loc[33]
    n.rubp = vars_gloss["New Name"].loc[34]
    n.ru5p = vars_gloss["New Name"].loc[35]
    n.pc_red = vars_gloss["New Name"].loc[36]
    n.fd_red = vars_gloss["New Name"].loc[37]
    n.nadp = vars_gloss["New Name"].loc[38]
    n.pi = vars_gloss["New Name"].loc[39]
    vars_gloss["New Name"].loc[40] # N_translocator does not exist
    n.mda = vars_gloss["New Name"].loc[41]
    n.h2o2 = vars_gloss["New Name"].loc[42]
    n.dha = vars_gloss["New Name"].loc[43]
    n.glutathion_ox = vars_gloss["New Name"].loc[44]
    n.tr_ox = vars_gloss["New Name"].loc[45]
    vars_gloss["New Name"].loc[46] #E_CBB_inactive
    n.ps2cs = vars_gloss["New Name"].loc[47]
    n.quencher = vars_gloss["New Name"].loc[48]
    vars_gloss["New Name"].loc[49] # PSI states
    n.lhcp = vars_gloss["New Name"].loc[50]
    n.tr_red = vars_gloss["New Name"].loc[51]
    vars_gloss["New Name"].loc[52] # E_CBB_active
    n.ascorbate = vars_gloss["New Name"].loc[53]
    n.glutathion_red = vars_gloss["New Name"].loc[54]
    n.a0 = vars_gloss["New Name"].loc[55]
    n.a1 = vars_gloss["New Name"].loc[56]
    n.a2 = vars_gloss["New Name"].loc[57]
    n.fluorescence = vars_gloss["New Name"].loc[58]
    vars_gloss["New Name"].loc[59]
    vars_gloss["New Name"].loc[60]
    vars_gloss["New Name"].loc[61]
    vars_gloss["New Name"].loc[62]
    vars_gloss["New Name"].loc[63]
    vars_gloss["New Name"].loc[64]
    vars_gloss["New Name"].loc[65]
    vars_gloss["New Name"].loc[66]
    vars_gloss["New Name"].loc[67]
    vars_gloss["New Name"].loc[68]
    vars_gloss["New Name"].loc[69]
    
    rates_gloss = pd.read_csv(rates_glossary_path, keep_default_na=False, dtype={"Glossary ID": "Int64"}, usecols=["Python Var", "Glossary ID"])
    
    rates_gloss = rates_gloss.set_index("Glossary ID")
    
    rates_gloss["New Name"] = rates_gloss["Python Var"].apply(new_name_creator)
    
    n.ps2 = rates_gloss["New Name"].loc[1]
    n.ptox = rates_gloss["New Name"].loc[2]
    n.atp_synthase = rates_gloss["New Name"].loc[3]
    rates_gloss["New Name"].loc[4] # ATPact
    n.proton_leak = rates_gloss["New Name"].loc[5]
    n.ex_atp = rates_gloss["New Name"].loc[6]
    rates_gloss["New Name"].loc[7] #v_Xcyc
    n.lhc_protonation = rates_gloss["New Name"].loc[8]
    n.b6f = rates_gloss["New Name"].loc[9]
    n.fnr = rates_gloss["New Name"].loc[10]
    rates_gloss["New Name"].loc[11] # v_FQR
    n.ndh = rates_gloss["New Name"].loc[12]
    n.cyclic_electron_flow = rates_gloss["New Name"].loc[13]
    n.lhc_state_transition_21 = rates_gloss["New Name"].loc[14]
    n.lhc_state_transition_12 = rates_gloss["New Name"].loc[15]
    n.violaxanthin_deepoxidase = rates_gloss["New Name"].loc[16]
    n.zeaxanthin_epoxidase = rates_gloss["New Name"].loc[17]
    n.lhc_deprotonation = rates_gloss["New Name"].loc[18]
    n.rubisco_carboxylase = rates_gloss["New Name"].loc[19]
    n.fbpase = rates_gloss["New Name"].loc[20]
    n.sbpase = rates_gloss["New Name"].loc[21]
    n.phosphoribulokinase = rates_gloss["New Name"].loc[22]
    n.ex_pga = rates_gloss["New Name"].loc[23]
    n.ex_dhap = rates_gloss["New Name"].loc[24]
    n.ex_gap = rates_gloss["New Name"].loc[25]
    n.ex_g1p = rates_gloss["New Name"].loc[26]
    n.phosphoglycerate_kinase = rates_gloss["New Name"].loc[27]
    n.gadph = rates_gloss["New Name"].loc[28]
    n.triose_phosphate_isomerase = rates_gloss["New Name"].loc[29]
    n.aldolase_dhap_gap = rates_gloss["New Name"].loc[30]
    n.transketolase_gap_f6p = rates_gloss["New Name"].loc[31]
    n.aldolase_dhap_e4p = rates_gloss["New Name"].loc[32]
    n.transketolase_gap_s7p = rates_gloss["New Name"].loc[33]
    n.ribose_phosphate_isomerase = rates_gloss["New Name"].loc[34]
    n.ribulose_phosphate_epimerase = rates_gloss["New Name"].loc[35]
    n.g6pi = rates_gloss["New Name"].loc[36]
    n.phosphoglucomutase = rates_gloss["New Name"].loc[37]
    n.ps1 = rates_gloss["New Name"].loc[38]
    n.ex_nadph = rates_gloss["New Name"].loc[39]
    n.ferredoxin_reductase = rates_gloss["New Name"].loc[40]
    n.ferredoxin_thioredoxin_reductase = rates_gloss["New Name"].loc[41]
    n.mda_reductase2 = rates_gloss["New Name"].loc[42]
    n.glutathion_reductase = rates_gloss["New Name"].loc[43]
    n.ascorbate_peroxidase = rates_gloss["New Name"].loc[44]
    n.mehler = rates_gloss["New Name"].loc[45]
    n.dehydroascorbate_reductase = rates_gloss["New Name"].loc[46]
    n.mda_reductase1 = rates_gloss["New Name"].loc[47]
    n.tr_activation = rates_gloss["New Name"].loc[48]
    n.tr_inactivation = rates_gloss["New Name"].loc[49]
    rates_gloss["New Name"].loc[50]
    rates_gloss["New Name"].loc[51]
    rates_gloss["New Name"].loc[52]
    rates_gloss["New Name"].loc[53]
    rates_gloss["New Name"].loc[54]
    rates_gloss["New Name"].loc[55]
    rates_gloss["New Name"].loc[56]
    rates_gloss["New Name"].loc[57]
    rates_gloss["New Name"].loc[58]
    rates_gloss["New Name"].loc[59]
    
    n.h = new_name_creator("H")
    n.co2 = new_name_creator("CO2")
    n.pi_ext = new_name_creator("Pi_ext")
    n.o2 = new_name_creator("O2")
    n.total_adenosines = new_name_creator("AP_tot")
    n.total_nadp = new_name_creator("NADP_tot")
    n.total_ferredoxin = new_name_creator("Fd_tot")
    n.e0 = new_name_creator("Enz0_")
    n.e = new_name_creator("E0_")
    n.total_ascorbate = new_name_creator("ASC_tot")

def uses_other_functions(func: Callable, called_funcs: None | dict = None, first_inspect_file: str | None = None) -> dict:
    if called_funcs is None:
        called_funcs = {}
        
    source = inspect.getsource(func)
    tree = ast.parse(source)
    func_globals = func.__globals__
    
    if first_inspect_file is None:
        first_inspect_file = inspect.getsourcefile(func)

    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            name = node.func.id
            if name in func_globals and callable(func_globals[name]):
                target_fn = func_globals[name]
                if inspect.getsourcefile(target_fn) == first_inspect_file and name not in called_funcs:
                    called_funcs[name] = target_fn
                    called_funcs = uses_other_functions(target_fn, called_funcs, first_inspect_file)
                    
    return called_funcs

def write_init(
    m: Model,
    model_name: str,
    path: Path = Path("__init__.py"),
) -> None:

    import_lines = f"""
from mxlpy import Model, Variable, InitialAssignment
from .derived_quantities import include_derived_quantities
from .rates import include_rates

__all__ = ["{model_name}"]\n
"""

    model_init = f"""
def {model_name}() -> Model:
    m = Model()

    m.add_parameters(
        {"{"}
"""
    
    for name, value in m._parameters.items():
        model_init += f'            "{name}": {value.value},\n'
        
    model_init += "        }\n    )\n\r    m.add_variables(\n        {\n"
    
    for name, value in m._variables.items():
        if isinstance(value.initial_value, InitialAssignment):
            other_funcs = uses_other_functions(value.initial_value.fn)
            for other_func_name, other_func in other_funcs.items():
                if f"def {other_func_name}(" not in import_lines:
                        import_lines += inspect.getsource(other_func) + "\n"
            import_lines += inspect.getsource(value.initial_value.fn) + "\n"
            val = f'InitialAssignment(fn={value.initial_value.fn.__name__}, args={value.initial_value.args}, unit="REPLACE")'
        else:
            val = f'Variable({value.initial_value}, unit="REPLACE")'
        model_init += f'            "{name}": {val},\n'
    
    model_init += "        }\n    )\n\r    m = include_derived_quantities(m)\n    m = include_rates(m)\n\n    return m"
    
    with Path.open(path, "w") as f:
        f.write(import_lines)
        f.write(model_init)



def write_reactions(
    m: Model,
    path: Path = "rates.py",
) -> None:
    import_lines = """
from mxlpy import Model, Derived
import numpy as np
from typing import cast
    """
    
    def_lines = ""
    
    add_reaction_msg = """
def include_rates(m: Model):
    """
    
    fns_set = set()
    
    passed_reacs = []
    
    for reac_name, reac in m._reactions.items():
        
        reac_fn_name = reac.fn.__name__
        
        if not hasattr(mxlbricks.fns, reac.fn.__name__):
            this_line = inspect.getsource(reac.fn) + "\n"
            if this_line not in def_lines:
                def_lines += this_line.replace(reac.fn.__name__, reac_fn_name)
        else:
            fns_set.add(reac.fn)
        
        stoichiometry = "{"
        for variable, stoic in reac.stoichiometry.items():
            if type(stoic) is Derived:
                if not hasattr(mxlbricks.fns, stoic.fn.__name__):
                    if f"def {stoic.fn.__name__}" not in def_lines:
                        def_lines += inspect.getsource(stoic.fn) + "\n"
                else:
                    fns_set.add(stoic.fn)
                stoic = f"Derived(fn={stoic.fn.__name__}, args={stoic.args}, unit={stoic.unit})"
                
            stoichiometry += f"\"{variable}\": {stoic}, "
        stoichiometry += "}"
        
        add_reaction_msg += f"""
    m.add_reaction(
        name=\"{reac_name}\",
        fn={reac_fn_name},
        args={reac.args},
        stoichiometry={stoichiometry}
    )"""
        
        for used_reac_name, used_reac in uses_other_functions(reac.fn).items():
             if hasattr(mxlbricks.fns, used_reac_name):
                 fns_set.add(used_reac)
        
    add_reaction_msg += "\n\n    return m"
    
    import_lines += """
from .basic_funcs import (
"""
    
    for fn in fns_set:
        import_lines += f"""    {fn.__name__},\n"""
    
    import_lines += ")\n\n"
        
    with open(path, "w") as f:
        f.write(import_lines)
        f.write(def_lines)
        f.write(add_reaction_msg)
        
    return

def write_derived(
    m: Model,
    path: Path = "derived_quantities.py",
) -> None:
    import_lines = """
from mxlpy import Model
from mxlpy.surrogates import qss
import numpy as np
import math
from typing import cast, Iterable
    """
    
    def_lines = ""
    
    add_derived_quantities_msg = """
def include_derived_quantities(m: Model):
    """
    
    fns_set = set()
    
    for derived_name, derived in m._derived.items():
        
        if not hasattr(mxlbricks.fns, derived.fn.__name__):
            other_funcs = uses_other_functions(derived.fn)
            if other_funcs != {}:
                
                for other_func_name, other_func in other_funcs.items():
                    if f"def {other_func_name}(" not in def_lines:
                        if not hasattr(mxlbricks.fns, other_func_name):
                            def_lines += inspect.getsource(other_func) + "\n"
                        else:
                            fns_set.add(other_func)
            def_lines += inspect.getsource(derived.fn) + "\n"
        else:
            fns_set.add(derived.fn)
        
        add_derived_quantities_msg += f"""
    m.add_derived(
        name=\"{derived_name}\",
        fn={derived.fn.__name__},
        args={derived.args},
    )
"""
    
    for surrogate_name, surrogate in m._surrogates.items():
        if not hasattr(mxlbricks.fns, surrogate.model.__name__):
            def_lines += inspect.getsource(surrogate.model) + "\n"
        else:
            fns_set.add(surrogate.fn)
            
        add_derived_quantities_msg += f"""
    m.add_surrogate(
        name="{surrogate_name}",
        surrogate=qss.Surrogate(
            model={surrogate.model.__name__},
            args={surrogate.args},
            outputs={surrogate.outputs}
        )
    )
"""
    for readout_name, readout in m._readouts.items():
        if not hasattr(mxlbricks.fns, readout.fn.__name__):
            def_lines += inspect.getsource(readout.fn) + "\n"
        else:
            fns_set.add(readout.fn)
            
        add_derived_quantities_msg += f"""
    m.add_readout(
        name="{readout_name}",
        fn={readout.fn.__name__},
        args={readout.args}
    )
"""
    add_derived_quantities_msg += "\n\n    return m"
    
    import_lines += """
from .basic_funcs import (
"""
    
    for fn in fns_set:
        import_lines += f"""    {fn.__name__},\n"""
    
    import_lines += ")\n\n"
        
    with open(path, "w") as f:
        f.write(import_lines)
        f.write(def_lines)
        f.write(add_derived_quantities_msg)
        
    return

def write_basic_funcs(
    m: Model,
    path: Path = Path("basic_funcs.py"),
) -> None:
    fns_set = set()

    for reac in m._reactions.values():
        if hasattr(mxlbricks.fns, reac.fn.__name__):
            fns_set.add(reac.fn)

        for stoic in reac.stoichiometry.values():
            if type(stoic) is Derived and hasattr(mxlbricks.fns, stoic.fn.__name__):
                fns_set.add(stoic.fn)
                
        for used_reac_name, used_reac in uses_other_functions(reac.fn).items():
             if hasattr(mxlbricks.fns, used_reac_name):
                 fns_set.add(used_reac)

    for derived in m._derived.values():
        if hasattr(mxlbricks.fns, derived.fn.__name__):
            fns_set.add(derived.fn)
            
    for readout in m._readouts.values():
        if hasattr(mxlbricks.fns, readout.fn.__name__):
            fns_set.add(readout.fn)

    msg = ""

    for fn in fns_set:
        msg += f"{inspect.getsource(fn)}\n"


    with Path.open(path, "w") as f:
        f.write(msg)
