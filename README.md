# EVolve
Linking planetary mantles to atmospheric chemistry through volcanism using EVo and FastChem.

# Overview #

EVolve is a linked mantle degassing and atmospheric growth code, which models the growth of a rocky planet's secondary atmosphere under the influence of volcanism.

# Installation #

EVolve is written in Python3, and is incompatible with Python 2.7. Two very useful tools to set up python environments:\
[Pip](https://pip.pypa.io/en/stable/) - package installer for Python\
[Anaconda](https://docs.continuum.io/) - virtual environment manager

1. Clone the repository **with submodules** and enter directory
   ```
   git clone --recurse-submodules git@github.com:pipliggins/EVolve.git
   ```
   **Note:** If you don't clone with submodules you won't get one of the modules used to run EVolve, the FastChem equilibrium chemistry code.
   
2. Compile [FastChem](https://github.com/exoclime/FastChem):
   ```
   cd fastchem
   git submodule update --init --recursive
   mkdir build
   cd build
   cmake -DUSE_PYTHON==ON ..
   make
   ```
   This will pull the pybind11 module required for the python bindings, and compile both the C++ code, and the python bindings which are used in EVolve to conect to FastChem.
   
   **Note:** FastChem is an external C++ module, used to compute atmospheric equilibrium chemistry. Therefore, to run on Windows, I recommend using WSL (Windows Subsystem for Linux) to make the process of compiling the C code easier. If you encounter installation issues relating to the cmake version, I found the [accepted answer here](https://askubuntu.com/questions/1203635/installing-latest-cmake-on-ubuntu-18-04-3-lts-run-via-wsl-openssl-error) to work for me. A list of the suggested terminal commands can also be found at the bottom of this README file.

   **Note:** FastChem currently cannot be compiled on Apple Macs which contain Apple silicon chips (M* processors) - https://github.com/NewStrangeWorlds/FastChem/issues/9 
  

3. Install dependencies using either Pip install or Anaconda. Check [requirements.txt](https://github.com/pipliggins/evolve/blob/master/requirements.txt) for full details.
   If using Pip, install all dependencies from the main directory of EVolve using
   ```
   pip3 install -r requirements.txt
   ```
    
# Running EVolve #

EVolve can be run either with or without the FastChem equilibrium chemistry in the atmosphere. To run Evolve with FastChem, from the main directory of EVolve run
```
python evolve.py inputs.yaml --fastchem
```
The available tags are:
* --fastchem   ).This will use fastchem to run equilibrium chemistry in the atmosphere, producing more chemical species than the magma degassing model uses and enabling the atmospheric equilibrium temperature to be lower than magmatic.

* --nocrust    ).This option stops a crustal reservoir from being formed out of the degassed melt which has been erupted. Instead, the degassed melt and any volatiles remaining in it are re-incorporated back into the mantle. If this tag is NOT used, the mantle mass will gradually reduce as there is no mechanism for re-introducing the crustal material back into the mantle implemented here.

All the input models for EVolve, and the submodules EVo and FastChem are stored in the 'inputs' folder: 

|     Filename    | Relevant module | Properties                                                                                                                                                                                                     |
|:---------------:|:---------------:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **atm.yaml**    | EVolve main     | Sets the pre-existing atmospheric chemistry and surface pressures + temperatures for the planet                                                                                                                |
| **mantle.yaml** | EVolve main     | Sets the initial planetary mantle/rocky body properties, including temperature, mass, fO2, the mantle volatile concentrations and the volcanic intrusive:extrusive ratio                                       |
| **planet.yaml** | EVolve main     | Sets generic planetary properties and important run settings, including planetary mass, radius, the amount of mantle melting occurring at each timestep and the size & number of timesteps the model will run. |
| chem.yaml       | EVo             | Contains the major oxide composition of the magma being input to EVo                                                                                                                                           |
| env.yaml        | EVo             | Contains the majority of the run settings and volatile contents for the EVo run.                                                                                                                               |
| output.yaml     | EVo             | Stops any graphical input from EVo compared to it's default settings                                                                                                                                           |

Files highlighted in bold should be edited by the user; all others are optimied for EVolve and/or will be edited by the code as it is running.
Explainations for each parameter setting in the EVolve files can be found at the bottom of this README file.

As EVolve runs, it creates and updates files in the outputs folder as follows:
| Filename           | Data                                                                                                                                                                         |
|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| atmosphere_out.csv | Planetary surface pressure and atmospheric composition for tracked molecules in units of volume mixing ratios (actually mo fraction), calculated after each time step                               |
| mantle_out.csv     | Mantle volatile budget and fO2 after each timestep                                                                                                                           |
| volc_out.csv       | The final pressure iteration from the EVo output file in each timestep (storing melt volatile contents, atomic volatile contents, gas speciation in mol & wt fractions, etc) |
| fc_input.csv         | Generated if fastchem is selected: The input to FastChem after atmospheric mixing, and hydrogen escape if that is occuring, for each timestep.    |
| fc_out.csv         | Generated if fastchem is selected: The results from FastChem after each timestep                                                             |



# Installation help for WSL #

If you see an error saying that the installed version of cmake is too low to install FastChem, try these commands:
**Please note this is just a suggestion based on what worked for me, try these workarounds at your own risk!**

```
sudo apt-get update
sudo apt-get install apt-transport-https ca-certificates gnupg software-properties-common wget

wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -

sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt-get update

sudo apt-get install cmake
```

