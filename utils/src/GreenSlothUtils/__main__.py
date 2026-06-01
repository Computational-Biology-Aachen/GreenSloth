from pathlib import Path, PurePath
from importlib import import_module
import importlib.util
import sys
from inspect import getmembers, isfunction

import click

from GreenSlothUtils import installerfuncs
from mxlpy import Model
from GreenSlothUtils import extract_select_to_gloss, check_gloss_to_model, update_from_main_gloss, write_python_from_gloss, export_glossselect_from_model, export_odes_as_latex


def test_for_dir(dire: Path):
    if not dire.is_dir():
        click.secho(f"The given path '{dire}' is not a working directory. Try Again!", fg="red", bold=True)
        dire = None
    
    return dire

def get_model(model_dir: Path, model_name: str) -> Model | None:
    spec = importlib.util.spec_from_file_location("model", model_dir / "model" / "__init__.py")
    if spec is None or spec.loader is None:
        msg = f"Cannot find spec for module {model_name}"
        raise ImportError(msg)

    module = importlib.util.module_from_spec(spec)
    sys.modules["model"] = module
    spec.loader.exec_module(module)

    for func in getmembers(module, isfunction):
        if func[0] == model_name:
            return func[1]()
    return None

def compare_gloss_to_model(model_dir: str | None, modelinfo_dir: str | None, modelgloss_dir: str | None) -> None:
    """Compares glosses to model"""
    # Check all the dirs
    model_dir = Path().absolute() if model_dir is None else Path(model_dir)
    
    model_dir = test_for_dir(model_dir)
    
    if model_dir is None:
        return
    
    modelinfo_dir = model_dir / "model_info" if modelinfo_dir is None else Path(modelinfo_dir)
    
    modelinfo_dir = test_for_dir(modelinfo_dir)
    
    if modelinfo_dir is None:
        return
    
    modelgloss_dir = modelinfo_dir / "model_glosses" if modelgloss_dir is None else Path(modelgloss_dir)
    
    modelgloss_dir = test_for_dir(modelgloss_dir)
    
    if modelgloss_dir is None:
        return
    
    click.secho("Checking for inconsistencies between model and glosses...", bg="white", fg="bright_black", bold=True)
        
    for i in ["comps", "rates", "params", "derived_comps", "derived_params"]:
        check_gloss_to_model(
            from_model=modelgloss_dir / f"model_{i}.csv",
            edit_gloss=modelinfo_dir / f"{i}.csv",
            check_col="Python Var",
            cli_flag=True
        )
        
#############################
# CLICK
#############################

@click.group()
def cli() -> None:
    pass

@cli.command()
@click.option("--path", "-p", default=None, help="Path to create model directory. Defaults to path here.")
@click.argument("model_name", metavar="<model-name>")
def initialize(path: str, model_name: str) -> None:
    """Create '<model-name>' directory."""
    try:
        installerfuncs.gs_install(model_name, path)
    except FileExistsError as exep:
        click.echo(str(exep))
        
    
@cli.command()
@click.option("--model-dir", "-md", default=None, help="Path to model directory. Defaults to path here")
@click.option("--modelinfo-dir", "-mid", default=None, help="Path to model info directory. Defaults to model-dir + 'model_info'")
@click.option("--modelgloss-dir", "-mgd", default=None, help="Path to where to store csvs. Defaults to model-dir + 'model_info/model_glosses/'")
@click.option("--extract-option", "-eo", default="all", help="""
              \b
              Parts of the model to extract. Possibilities:
              \t 'all',
              \t 'variables',
              \t 'parameters',
              \t 'derived_variables',
              \t 'derived_parameters',
              \t 'reactions',
              \t [default: 'all']
            """)
@click.option("--check/--no-check", default=True, help="Check for inconsistencies with 'compare_gloss_to_model'")
def from_model_to_gloss(model_dir: str | None, modelinfo_dir: str | None, modelgloss_dir: str | None, extract_option: str, check: bool) -> None:
    """Generate temporary Glosses from model info."""
    
    # Check all the dirs
    model_dir = Path().absolute() if model_dir is None else Path(model_dir)
    
    model_dir = test_for_dir(model_dir)
    
    if model_dir is None:
        return
    
    modelinfo_dir = model_dir / "model_info" if modelinfo_dir is None else Path(modelinfo_dir)
    
    modelinfo_dir = test_for_dir(modelinfo_dir)
    
    if modelinfo_dir is None:
        return
    
    modelgloss_dir = modelinfo_dir / "model_glosses" if modelgloss_dir is None else Path(modelgloss_dir)
    
    modelgloss_dir = test_for_dir(modelgloss_dir)
    
    if modelgloss_dir is None:
        return
    
    model_name = model_dir.name
    
    if extract_option.lower() not in ["all", "variables", "parameters", "derived_variables", "derived_parameters", "reactions"]:
        click.echo(f"The given extraction option '{extract_option}' is not valid!")
        return
    
    try:
        m = get_model(model_dir, model_name)
    except ImportError:
        click.secho("Could not import Model! Are the directory and name correct?", fg="red", bold=True)
        return
    
    if extract_option.lower() in ["all", "variables"]:
        extract_select_to_gloss(
            select=m.get_raw_variables(),
            column_names=[
                "Name",
                "Common Abbr.",
                "Paper Abbr.",
                "KEGG ID",
                "Python Var",
                "Glossary ID",
            ],
            pythonvar_col="Python Var",
            path_to_write= modelgloss_dir / "model_comps.txt",
        )

    if extract_option.lower() in ["all", "parameters"]:
        extract_select_to_gloss(
            select=m.get_parameter_values(),
            column_names=[
                "Short Description",
                "Common Abbr.",
                "Paper Abbr.",
                "Value",
                "Unit",
                "Python Var",
                "Reference",
            ],
            pythonvar_col="Python Var",
            path_to_write= modelgloss_dir / "model_params.txt",
            value_col="Value"
        )
        
    if extract_option.lower() in ["all", "derived_variables"]:
        surr_outputs_dict = {}
        for surr in m.get_raw_surrogates().values():
            for output in surr.outputs:
                surr_outputs_dict[output] = None
        full_dict = m.get_derived_variables() | surr_outputs_dict | m.get_raw_readouts()
        extract_select_to_gloss(
            select=full_dict,
            column_names=[
                "Name",
                "Common Abbr.",
                "Paper Abbr.",
                "KEGG ID",
                "Python Var",
                "Glossary ID",
            ],
            pythonvar_col="Python Var",
            path_to_write= modelgloss_dir / "model_derived_comps.txt",
        )
        
    if extract_option.lower() in ["all", "derived_parameters"]:
        extract_select_to_gloss(
            select=m.get_derived_parameters(),
            column_names=["Short Description", "Common Abbr.", "Paper Abbr.", "Python Var"],
            pythonvar_col="Python Var",
            path_to_write= modelgloss_dir / "model_derived_params.txt",
        )
        
    if extract_option.lower() in ["all", "reactions"]:
        extract_select_to_gloss(
            select=m._reactions,
            column_names=[
                "Short Description",
                "Common Abbr.",
                "Paper Abbr.",
                "KEGG ID",
                "Python Var",
                "Glossary ID",
            ],
            pythonvar_col="Python Var",
            path_to_write= modelgloss_dir / "model_rates.txt",
        )
        
    if check:
        compare_gloss_to_model(model_dir, modelinfo_dir, modelgloss_dir)

@cli.command()
@click.option("--model-dir", "-md", default=None, help="Path to model directory. Defaults to path here")
@click.option("--modelinfo-dir", "-mid", default=None, help="Path to model info directory. Defaults to model-dir + 'model_info'")
@click.option("--modelgloss-dir", "-mgd", default=None, help="Path to model glosses are stored. Defaults to model-dir + 'model_info/model_to_glosses/'")
def compare_gloss_to_model_command(model_dir: str | None, modelinfo_dir: str | None, modelgloss_dir: str | None) -> None:
    """Compares glosses to model"""
    compare_gloss_to_model(model_dir, modelinfo_dir, modelgloss_dir)
    
@cli.command()
@click.option("--maingloss-dir", "-magd", default=None, help="Path to directory with main gloss. Defaults to parent of here.")
@click.option("--model-dir", "-md", default=None, help="Path to model directory. Defaults to path here")
@click.option("--modelinfo-dir", "-mid", default=None, help="Path to model info directory. Defaults to model-dir + 'model_info'")
@click.option("--add/--no-add", default=False, help="Add new entries to main gloss. Defaults to False.")
def update_glosses_from_main(maingloss_dir: str | None, model_dir: str | None, modelinfo_dir: str | None, add: bool) -> None:
    """Update glosses from main"""
    
    # Check all the dirs
    maingloss_dir = Path().absolute().parent if maingloss_dir is None else Path(maingloss_dir)
    
    maingloss_dir = test_for_dir(maingloss_dir)
    
    if maingloss_dir is None:
        return
    
    model_dir = Path().absolute() if model_dir is None else Path(model_dir)
    
    model_dir = test_for_dir(model_dir)
    
    if model_dir is None:
        return
    
    modelinfo_dir = model_dir / "model_info" if modelinfo_dir is None else Path(modelinfo_dir)
    
    modelinfo_dir = test_for_dir(modelinfo_dir)
    
    if modelinfo_dir is None:
        return
    
    model_name = model_dir.name
    
    for i, j in zip(
        ("comps", "comps", "rates"), ("comps", "derived_comps", "rates"), strict=False
    ):
        update_from_main_gloss(
            main_gloss_path=maingloss_dir / f"{i}_glossary.csv",
            gloss_path=modelinfo_dir / f"{j}.csv",
            add_to_main=add,
            model_title=model_name,
            cli_flag=True
        )
        
@cli.command()
@click.option("--model-dir", "-md", default=None, help="Path to model directory. Defaults to path here")
@click.option("--modelinfo-dir", "-mid", default=None, help="Path to model info directory. Defaults to -md + 'model_info'")
@click.option("--glosstopython-dir", "-gpd", default=None, help="Path to gloss to python directory. Defaults to -mid + 'python_written/gloss_to_python'")
def python_from_gloss(model_dir: str | None, modelinfo_dir: str | None, glosstopython_dir: str | None) -> None:
    """Write Python Variables from Glossaries"""
    
    # Check all the dirs
    model_dir = Path().absolute() if model_dir is None else Path(model_dir)
    
    model_dir = test_for_dir(model_dir)
    
    if model_dir is None:
        return
    
    modelinfo_dir = model_dir / "model_info" if modelinfo_dir is None else Path(modelinfo_dir)
    
    modelinfo_dir = test_for_dir(modelinfo_dir)
    
    if modelinfo_dir is None:
        return
    
    glosstopython_dir = modelinfo_dir / "python_written/gloss_to_python" if glosstopython_dir is None else Path(glosstopython_dir)
    
    glosstopython_dir = test_for_dir(glosstopython_dir)
    
    if glosstopython_dir is None:
        return
    for i in ['comps', 'rates', 'params', 'derived_comps', 'derived_params']:
        write_python_from_gloss(
            path_to_write=glosstopython_dir / f'{i}.txt',
            path_to_glass=modelinfo_dir / f'{i}.csv',
            var_list_name=f'{i}_table'
        )
        
@cli.command()
@click.option("--model-dir", "-md", default=None, help="Path to model directory. Defaults to path here")
@click.option("--modelinfo-dir", "-mid", default=None, help="Path to model info directory. Defaults to -md + 'model_info'")
@click.option("--modeltolatex-dir", "-mld", default=None, help="Path to model to latex directory. Defaults to -mid + 'python_written/model_to_latex'")
def latex_from_model(model_dir: str | None, modelinfo_dir: str | None, modeltolatex_dir: str | None) -> None:
    """Write LaTex from Model"""
    
    # Check all the dirs
    model_dir = Path().absolute() if model_dir is None else Path(model_dir)
    
    model_dir = test_for_dir(model_dir)
    
    if model_dir is None:
        return
    
    modelinfo_dir = model_dir / "model_info" if modelinfo_dir is None else Path(modelinfo_dir)
    
    modelinfo_dir = test_for_dir(modelinfo_dir)
    
    if modelinfo_dir is None:
        return
    
    modeltolatex_dir = modelinfo_dir / "python_written/model_to_latex" if modeltolatex_dir is None else Path(modeltolatex_dir)
    
    modeltolatex_dir = test_for_dir(modeltolatex_dir)
    
    if modeltolatex_dir is None:
        return
    
    model_name = model_dir.name
    
    try:
        m = get_model(model_dir, model_name)
    except ImportError:
        click.secho("Could not import Model! Are the directory and name correct?", fg="red", bold=True)
        return
    
    for i in ["rates", "derived_comps", "derived_params"]:
        export_glossselect_from_model(
            m=m,
            write_path=modeltolatex_dir / f"{i}.txt",
            gloss_path=modelinfo_dir / f"{i}.csv",
        )
        
    export_odes_as_latex(m=m, path_to_write=modeltolatex_dir / "model_odes.txt")
    
def start_cli() -> None:
    cli(max_content_width=200)