# Mie solver implementation in python
This is the Mie solver implementation using the miepython package

## Installation
To run locally just download the github repo, open the terminal in the project folder and depending on what you use for package management, create an environment with the dependencies and run the files.

### Using Pip
Open the terminal at the project folder location and create a new environment:

`python -m venv <environment_name>`

Then activate the environment and install the packages running the following commang:

`pip install -r requirements.txt`

### Using conda
Open the conda terminal in the project folder and build the environment using the yml file:

`conda env create -f environment.yml`

After that you can activate the environment (it is called "scat"):

`conda activate scat`

You should be able to see this environment on your IDE to compile the files and can select the scat environment kernel inside the jupyter notebooks.

## Project Structure
The project has the following structure:
- The main.py file has a simple script that shows how to use the miepython package, as a sample I used gold nanoparticles and it seems to work without issues.
- The test_miepython_interpolation.ipynb is a jupyter notebook where I tested some of my interpolation functions to see if they would work properly and would not screw up the calculation, they seem to be working without issues.
- The test_silicon_for_comparison.ipynb is a jupyter notebook where I try to replicate the CSV file you sent me with the different efficiencies, I used the absorption one as a base, this notebook at the will generate a similar CSV but with the miepython calculation.
- Inside the tests folder, there are both the CSVs from the mathematica implementation and the miepython implementation, and a jupyter notebook where it loads both CSVs and plots the efficiencies for the same radii.

## Results so far
Right now miepython seems to be working properly, some of the efficiencies do not look that different like the one for 36 nm 

<img width="565" height="480" alt="Abs_36nm" src="https://github.com/user-attachments/assets/c723a9b2-89a2-4019-a329-290c746b03eb" />


it sort of follows this trend until the 70 nm marks where the differences start to show

<img width="570" height="480" alt="Abs_71nm" src="https://github.com/user-attachments/assets/3c18e087-de39-4611-a578-99e7890e1959" />

Though there are some sizes that still have some similarities, you can generate these plots and play around them in the check.ipynb notebook, and I've also added some of the graphs for a few sizes in the results folder.

# Installing requirments

To run `mie_solver_collab`, you need to set up a dedicated Conda environment. The following steps will create the environment with Python 3.11 and install all required packages using the `environment.yml` file.

1. Create the Conda environment and install the dependencies:
   ```bash
   conda env create -f environment.yml -n conda_mie_solver
   ```

2. Activate the new environment:
   ```bash
   conda activate conda_mie_solver
   ```