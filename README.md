# GreenSloth

**A curated, open repository of dynamic photosynthesis models.**

---

## Why GreenSloth?

Photosynthesis underpins nearly all life on Earth, yet the models that describe it are scattered across supplementary materials, personal repositories, and institutional pages that disappear when funding ends. Despite generous funding and overwhelming global effort, the field has no shared infrastructure in terms of models. This leads to similar models being redeeloped, stalling collaborative progress.

GreenSloth was built to change that.

Inspired by the legacy of E-photosynthesis and the growing gap between the pace of model development and the community's ability to access, compare, and build on that work, GreenSloth is a living repository of curated, documented, and modular photosynthesis models. From the classical steady-state Farquhar–von Caemmerer–Berry (FvCB) framework to state-of-the-art kinetic and dynamic systems models, including e-photosynthesis.

---

## Who is GreenSloth for?

**Experimentalists**  
You collected beautiful data — chlorophyll fluorescence transients, 
gas exchange curves, proteomics under fluctuating light — and you want 
to know whether the current mechanistic understanding of photosynthesis 
can reproduce what you see. GreenSloth gives you a curated entry point 
into simulation tools, with enough documentation to connect your 
experimental conditions to the right model.

**Computational Biologists & Modellers**  
Building a model of a photosynthetic subsystem shouldn't mean starting 
from scratch. GreenSloth organises models into modular, interoperable 
components — electron transport, carbon fixation, photoprotection, 
stomatal regulation — so you can assemble, compare, and extend existing 
work rather than duplicating it.

**AI & Machine Learning Researchers**  
Hybrid approaches combining data-driven methods with first-principles 
models are among the most promising frontiers in systems biology. 
GreenSloth provides the mechanistic backbone: curated models that encode 
decades of domain knowledge, ready to serve as physics-informed priors, 
training constraints, or ground-truth benchmarks for neural network 
architectures.

---

## What's in the repository?

- Curated implementations of landmark photosynthesis models
- Standardised metadata: biological scope, timescale, species, 
  experimental context, and known limitations
- Modular components tagged by subsystem
- Links to primary literature and original datasets where available
- Roadmap for community contributions

---

## The name

The sloth is the only known mammal whose fur hosts photosynthetic organisms. Slow, deliberate, and surprisingly green, a fitting mascot for out project.

---
## What is this

All the models in the ecosystem of GreenSloth are found in this repository. Each model has its own directory, with a specific structure. The overarching name is the last Name of the first author and the date of publication of the model. This directory is best created using the [GreenSlothUtils](https://github.com/ElouenCorvest/GreenSlothUtils), as it will automatically use the name provided and insert it in the right places. Inside the model directory, you can find the following:

```bash
Corvest2000
├── figures
│   ├── demonstrations.ipynb
│   └── paper_figures.ipynb
├── model  
│   ├── __init__.py
│   ├── derived_quantities.py
│   ├── rates.py
│   └── basic_funcs.py
├── model_info
│   ├── comps.csv
│   ├── derived_comps.csv
│   ├── derived_params.csv
│   ├── params.csv
│   ├── rates.csv
│   ├── model_glosses
│   └── python_written
│       ├── gloss_to_python
│       └── model_to_latex
│ 
└── README_script.py
```

The `figures` directory includes the template Jupyter Notebooks for the recreations of the figures of the original paper and for demonstrations of the model.

The `model` directory includes the Python code of the model, which is seperated into different files for better readability. The `__init__.py` file is the main file, which imports the other files and contains the main function of the model. The `derived_quantities.py` file contains the derived quantities of the model, and the `rates.py` file contains the rates of the model. The `basic_funcs.py` file contains the basic functions of the model, which are used in the other files.

The `model_info` directory contains all the information about the model. The `.csv` files need to be filled with the information of the model, which is used in the README file. The comps and rates can be helped by the addition of the Glossary IDs. The `model_glosses` directory is filled after usage of the `GreenSlothUtils` to create the glossaries extracted from the model. The `python_written` directory contains pointers to be copied over to the README script, which are created by `GreenSlothUtils` from the model. 

The `README_script.py` file is the script that generates the README file of the model, which is the main file that is shown on the website. This script needs to be filled with the python variables, latex equations, summary, figure recreations, and demontrations of the model. The script is then run, which generates the README file.

## How to Contribute

GreenSloth started as a master thesis project of M.Sc. Elouen Corvest who designed its first fully functioning prototype. After his defense, the project was expanded by the Matuszynska Lab and now grows through the community. If you have a model, a correction, or a missing reference, see CONTRIBUTING.md.

To create a new model directory to be then included here, the best way is to follow the instructions in the [GreenSlothUtils](https://github.com/ElouenCorvest/GreenSlothUtils). This will automatically create the directory with the right structure and also fill in some of the files with the right information. After that, you can fill in the rest of the files with the information of the model. The README script is the most important file, as it generates the README file that is shown on the website. Therefore, it needs to be filled with care and attention to detail.

Once the model directory is created and everything is filled out, you can create a pull request to this repository, which whill then be reviewed and merged by the maintainers of this repository. If your model directory follows the same format as the other model directories, it can be accepted. if it does not follow the same format, it will be sent back to you with comments on what needs to be changed.

After acceptance, the maintainers will then add the model to the website, which can take some time.
