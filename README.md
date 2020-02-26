# pyOrigamiEM
DNA Origami Goniometer assisted cryoEM data analysis package. Python scripts and classes in this package can be used for other cryo-EM projects. Please cite ... if you use this python package for your work. 

Part of cryosparc `.cs` to relion `.star`data conversion script is adapted from Daniel Asarnow's [pyem python package](https://github.com/asarnow/pyem).

## Dependencies

- [Relion](https://github.com/3dem/relion)
- Python 3.5 or above. 

## Installation

Installing this package inside python virtual environment is high encouraged. After installing `virtualenv` and `virtualenvwrapper`, create a python3 virtual environnment.

`>mkvirtualenv em -p python3`.

Remember to activate the virtualenvironment

`>workon em`

For most up to date version of the package, clone or download. Inside the `cryoorigami` package folder, execute

`> (em) pip install .` 

If you would like to run the scripts in any path, make sure to add the `bin` folder in `$PATH` variable in your `.bash_profile` or `.bashrc` file.

## Descriptions for commonly used scripts

`em_addstarcols.py`: Add new columns to a star file.

`em_align2D.py`: Align particles from 2D classification output to a reference class average.

`em_alignclassaverages.py`: Align class averages to a single class.   


`
...
