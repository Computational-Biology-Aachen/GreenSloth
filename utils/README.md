# GreenSlothUtils

This package includes all the commands needed to simplify the work on the GreenSloth Project.

## How to Create a model for GreenSloth (Using CLI)

In thsi section you are getting a step-by-step guide on how the easiest way is to prepare a model for GreenSloth.

## Installation



### 0. All the commands

Using `GreenSlothUtils` enables an easier way to contribute to GreenSloth. While there are many different functions and upsides this custom package gives, the most important ones are summarized in a CLI interface. These commands initializes your work for contributing a model to GreenSloth.

Once you have the `GreenSlothUtils` package in your environment (we recommend using [`pixi`](https://pixi.sh/latest/)) you can enter the following in your Terminal. Any Terminal should work.

```bash
GreenSloth-init --help
```

This command shows you all the availeble options in the CLI interface. All the commands are also explained here and given an order to fully encapsulate the way, we best believe, you should summarise for GreenSloth.

### 1. Create the Directory

Every model should in a seperate directory, for ease of distinction. by using the following command you can do just that:

```bash
GreenSloth-init initialize <model-name> -p=<path>
```

The argument `<model-name>` is mandatory and represents the name of the directory. We recommend using the name of your model, as it makes it much easier. Additionally, all created files that use the model name will also be using this given argument. So choose wisely!

The option `-p=<path>`, is optional and lets you give the command a specific path where the directory should be created. If this is ommited, the path where the Terminal points to will be taken.

### 2. Create your model using `MxlPy`

This step is sadly not easily replaced by `GreenSlothUtils`. You will have to create your own working model and we recommend doing so with `MxlPy` if applicable. To see how to use that package please go its [docs](https://github.com/Computational-Biology-Aachen/MxlPy). Inside your newly created model directory, you will see the model directory, that contains several python-scripts. These scripts are there to guide you through the model-creating process and help you to keep a tidy representation fo your model.

The `__init__.py` file is the main file of your model. It collects all the other files and actually combines them to create your model, which should be using the name of the parent directory. The only things you should add to this file, are the variables and the parameters.

The `basic_funcs.py` file is empty, as it is a file for you to put common functions used by your model. The rest of the model files shoudl automatically import all the functions you define here.

The `derived_quantities.py` and the `rates.py` file should include, as their appropriate names show, the derived quantitites and the rates of your model. Please include any addition to these files inside the function, which name's you should not change! These two functions are directly called by the `__init__.py` file to construct the rest of your model.

At the end, you should be able to import the model folder and use the `<model_name>()` function (where you replace the `<model_name>` with the appropriate name) to constuct your entire model in any other python-script. This function is also in some of the commands used by `GreenSlothUtils` 

### 3. Extract Infromation from your Model

To fascilitate the summary of your model, `GreenSlothUtils` includes a function that extracts common model information, like the variables, parameters, both their derived parts, and the rates. This information is inputted into seperated csv tables that follow the same format as the overall glossaries.

```bash
Usage: GreenSloth-init from-model-to-gloss [OPTIONS]

  Generate temporary Glosses from model info.

Options:
  -md, --model-dir TEXT        Path to model directory. Defaults to path here
  -mid, --modelinfo-dir TEXT   Path to model info directory. Defaults to
                               model-dir + 'model_info'
  -mgd, --modelgloss-dir TEXT  Path to where to store csvs. Defaults to model-
                               dir + 'model_info/model_to_glosses/'
  -eo, --extract-option TEXT   Parts of the model to extract. Possibilities:
                                  'all',
                                  'variables',
                                  'parameters',
                                  'derived_variables',
                                  'derived_parameters',
                                  'reactions',
                                  [default: 'all']
  --check / --no-check         Check for inconsistencies with
                               'compare_gloss_to_model'
```

`-md, --model-dir`

The path to the model directory, usually created by `GreenSloth-init initialize`. Defaults to path of execution.

`-mid, --modelinfo-dir`

The path to the directory where you wish to store all the model information. Defaults to `-md` + `'model_info'`.

`-mgd, --modelgloss-dir`

The path where to store the created csv files. Defaults to `-mid` + `'model_to_glosses'` 

`-eo, --extract-options`

A string of what should be extracted from the model. The possible entries are: `all`, `variables`, `parameters`, `derived_variables`, `derived_parameters`, and `reactions`. Defaults to `all`.

`--check / --no-check`

Boolean to compare the created csvs with the general glossaries, using the function [`compare-gloss-to-model`](). Defaults to `--check`.

### 4. Correct and Use Model Info

Once you have extracted all the model information, you should now fill in the main glossaries created by `GreenSloth-init initialize` appropriately. While doing so, check if some information is also contained in the main overarching glossaries, if so you can just input the Glossary ID and the abbreviation used in the paper, and later code will fill it in for you. 

### 5. Update from Main Glossaries

Once you have sieved through your entire model information to and compiled your glossaries, you may want to update them to fit the main overarching glossaries. This will look at the IDs you have given, and will all the rest of the information. Optionally, you can also update the main glossaries with new entries that you did not give an ID to.

```bash
Usage: GreenSloth-init update-glosses-from-main [OPTIONS]

  Update glosses from main

Options:
  -magd, --maingloss-dir TEXT  Path to directory with main gloss. Defaults to
                               parent of here.
  -md, --model-dir TEXT        Path to model directory. Defaults to path here
  -mid, --modelinfo-dir TEXT   Path to model info directory. Defaults to
                               model-dir + 'model_info'
  --add / --no-add             Add new entries to main gloss. Defaults to
                               False.
```

`-magd, --maingloss-dir`

Path to the directory containing the main overarching glossaries. Defaults to parent of path of execution.

`-md, --model-dir`

The path to the model directory, usually created by `GreenSloth-init initialize`. Defaults to path of execution.

`-mid, --modelinfo-dir`

The path to the directory where you wish to store all the model information. Defaults to `-md` + `'model_info'`.

`--add / --no-add`

Boolean value to add new entries to the main overarching glossary. Defaults to `--no-add`. We suggest to do this step at the end and only if you are absolutely sure, the new entries are actually new, as correcting the main glossary can be a pain.

### 6. Correct your model

Once you have updated your glossaries, you may need to change the variables used in your python model. To help you with that, you can use the `--check` option from `GreenSloth-init from-model-to-gloss` that will tell you if there are any inconsistencies. Keep in mind, that the different categories of information may not reflect your way of seperating them completely. For example, a model may have protons as a variable, while another has it as a parameter. Therefore calculating the pH could result in a derived variable or parameter. Luckily, the check only gives a you a list of inconsistencies, therefore it should be very easy to find the culrpits and determine how and if they need to be corrected.

### 7. Convert to Python variables

As the summary of your model should be written in a markdown file, we need to write the markdown in a Python script. This allows us to use variables and functions to stay consistent in the entire file. One big adavantage is that we can store the math expression of all our model variables, parameters, and reactions inside seperate python variables. And `GreenSlothUtils` enables that automatically.

```bash
Usage: GreenSloth-init python-from-gloss [OPTIONS]

  Write Python Variables from Glossaries

Options:
  -md, --model-dir TEXT           Path to model directory. Defaults to path here
  -mid, --modelinfo-dir TEXT      Path to model info directory. Defaults to -md + 'model_info'
  -gpd, --glosstopython-dir TEXT  Path to gloss to python directory. Defaults to -mid + 'python_written/gloss_to_python'
```

Several .txt files will be created that include the appropriate stored python variables, that can be copied over to the `README_script.py`. It checks if the .txt files that are already there, need to be changed, and if so will generate an update section. So you can see when something was changed in these files.

`-md, --model-dir`

The path to the model directory, usually created by `GreenSloth-init initialize`. Defaults to path of execution.

`-mid, --modelinfo-dir`

The path to the directory where you wish to store all the model information. Defaults to `-md` + `'model_info'`.

`-gpd, --glosstopython-dir`

The path to the directory you wish to store the python in a .txt file. Defaults to `-mid` + `'python_written/gloss_to_python'`.

### 8. Create Latex from Model

In the summary of your model, you will also need to add the LaTex version of your reactions, derivedd equations, and ODE system. With `GreenSlothUtils` there is an easy way to do so. `mxlPy` has its own intern way to give you acces to LaTex version of your model, but as we are working with the prior created python variables and want some custom formatting, `GreenSlothUtils` uses its own methods.

```bash
Usage: GreenSloth-init latex-from-model [OPTIONS]

  Write LaTex from Model

Options:
  -md, --model-dir TEXT          Path to model directory. Defaults to path here
  -mid, --modelinfo-dir TEXT     Path to model info directory. Defaults to -md + 'model_info'
  -mld, --modeltolatex-dir TEXT  Path to model to latex directory. Defaults to -mid + 'python_written/model_to_latex'
```

Several .txt files will be created that include the appropriate stored LaTex version of your model using the python variables, that can be copied over to the `README_script.py`. It checks if the .txt files that are already there, need to be changed, and if so will generate an update section. So you can see when something was changed in these files.

`-md, --model-dir`

The path to the model directory, usually created by `GreenSloth-init initialize`. Defaults to path of execution.

`-mid, --modelinfo-dir`

The path to the directory where you wish to store all the model information. Defaults to `-md` + `'model_info'`.

`-mld, --modeltolatex-dir`

The path to the directory you wish to store the LaTex in a .txt file. Defaults to `-mid` + `'python_written/gloss_to_python'`.

### 9. Finishing touches

Once you have copied everything into your `README_script.py` you can begin filling out the other information, like a brief a summary and the recreation of the figures of the paper. Once these things are done, you have finished with summarizing your model, and it can be inserted into the website! 