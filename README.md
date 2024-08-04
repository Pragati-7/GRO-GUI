# GRO-GUI: An Interactive Graphical User Interface for GROMACS

## GROMACS Lysozyme in Water Tutorial
This repository provides a tutorial for running a molecular dynamics simulation of lysozyme in water using GROMACS, a versatile package for performing molecular dynamics.

## Table of Contents

1. [Introduction](#introduction)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Usuage](#usuage)
5. [Workflow Steps](#workflow-steps)
6. [Trobleshooting](troubleshooting)
7. [Contributing](#contributing)

## Introduction


## Requirements

* Operating System: Linux
* GROMACS: Ensure GROMACS is installed and configured in your environment. You can follow the official [GROMACS installation guide](https://manual.gromacs.org/documentation/current/install-guide/index.html).
* Python: Ensure Python is installed (preferably version 3.6 or above).
* Dependencies: The following Python packages are required:
    * tkinter
    * Pillow

You can install the Python dependencies using:

```bash
pip install Pillow
pip install tkinter
```

## Installation

1. Clone the repository:

   To clone the given repository-

   ```
   git clone https://github.com/Pragati-7/GRO-GUI.git
   cd GRO-GUI
   ```

3. Set up GROMACS:

   Ensure GROMACS is properly installed and added to your system path. You can check your installation by running:

   ```
   gmx --version
   ```

## Usuage

1. Run the GUI:

```
python main.py
```

2. Follow the Workflow:

Use the GUI to proceed through each step of the simulation, as outlined below.

## Workflow Steps

1. Generate Topology:

   * Script: 1Generate_topology.py
   * Description: Generates the topology files necessary for simulation.

2. Define Box and Solvate:

   * Script: 2Define_box_and_solvate.py
   * Description: Defines the simulation box and adds water molecules.

3. Add Ions:

   * Script: 3Adding_ions.py
   * Description: Neutralizes the system by adding ions.

4. Energy Minimization:

   * Script: 4Energy_minimization.py
   * Description: Minimizes the potential energy of the system.
     
5. Equilibration:

   * Script: 5Equilibration.py
   * Description: Equilibrates the system to prepare for production MD.
   
6. Production MD:

    * Script: 6Production_MD.py
    * Description: Performs the production molecular dynamics simulation.

7. Analysis:

    * Script: 7Analysis.py
    * Description: Analyzes the results of the simulation.

## Troubleshooting

* GROMACS not found: Ensure GROMACS is installed and the gmx command is in your system's PATH.
* Python errors: Check that all required Python packages are installed.

## Contributings

For any issues or inquiries, please contact [Pragati Singh] at [pragatisinghsingh45@gmail.com].
