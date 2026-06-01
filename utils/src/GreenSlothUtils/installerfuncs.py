import os
from pathlib import Path
from typing import Optional
import click

def iterate_files(path_to_scan: Path, model_name: str, target_dir: Path, dirs = "") -> None:

    for src_file in os.scandir(path_to_scan):
        if os.path.isdir(src_file):
            new_dirs = dirs
            new_dirs += src_file.name + '/'

            os.makedirs(
                target_dir / new_dirs,
                exist_ok=True
            )

            iterate_files(src_file, model_name, target_dir, dirs=new_dirs)

        elif os.path.isfile(src_file):
            with open(src_file, 'r') as f:
                content = f.read()
            if src_file.name.endswith(('.py', ".ipynb")):
                content = content.replace('{{MODEL_NAME}}', model_name)

            dest = target_dir / dirs / src_file.name

            with open(dest, 'w') as f:
                f.write(content)


def gs_install(
    model_name: str,
    target_dir: Optional[Path] = None,
):
    if target_dir is None:
        target_dir = Path('./')

    try:
        os.makedirs(
            target_dir / model_name,
            exist_ok=False
        )
    except FileExistsError:
        raise FileExistsError(f"There already exists a folder with name '{model_name}'")

    iterate_files(Path(__file__).parent / 'modelinit', model_name, target_dir= target_dir / model_name)
    
    for i in ["model_glosses", "python_written/gloss_to_python", "python_written/model_to_latex"]:
        try:
            Path.mkdir(target_dir / model_name / "model_info" / i, parents=True, exist_ok=True)
        except FileExistsError:
            raise FileExistsError(f"Something went wrong!")