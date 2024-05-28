
# Interdigital Capacitor (IDC) Capacitance Calculation, Simulation and Testing. 
This repository is part of the laboratory project for the Laboratory of Sensors and Sensing Systems course in the Master's degree program in 
Automation Engineering and Control of Complex Systems. It contains MATLAB scripts for calculating and simulating the capacitance of a 
3-layer interdigital capacitor (IDC) using various models and methods. Additionally, it includes datasets and analyses of 
analog measurements conducted in the lab.
The files included in this repository are:

1. `c_idc3k.m`
2. `idcSimulation.m`
3. `bareSensorMeasures.m`
4. `electronicMeasures.m`

## Files Overview

### 1. `c_idc3k.m`

This script contains the function `c_idc3k` which calculates the capacitance of a 3-layer IDC using the Kim model. 

**Input Parameters:**
- `eps1`, `eps2`, `eps3`: Dielectric permittivity of layers 1, 2, and 3, respectively.
- `h1`, `h2`, `h3`: Thickness of layers 1, 2, and 3, respectively.
- `b`: Finger width.
- `d`: Finger spacing.
- `l`: Overlapping finger length.
- `n`: Number of IDC finger pairs.
- `display` (optional): Boolean to display intermediate capacitance values.

**Output Parameter:**
- `Cidc`: Overall IDC Capacitance in Farads (F).

### 2. `idcSimulation.m`

This script defines the parameters for simulating the Interdigital Capacitor (IDC) and compares the results with the interdigital model provided by the RF and Mixed Signal MATLAB module.
**Features:**
- Calculation of IDC capacitance using the Kim model.
- Comparison with the interdigital MATLAB model.
- Sensitivity analysis of the capacitance with respect to variations in layer thickness (`h2`) and permittivity (`eps2`).
- Plotting the results for visualization.

### 3. `bareSensorMeasures.m`

This script analyzes the resistance and capacitance measurements of various sensor configurations. 

**Features:**
- Bar plots for resistance and capacitance (mean and standard deviation).
- Analysis of sensor measures using inkjet printing camera data.
- Comparison of capacitance values obtained with set and printed parameters.
- Finetuning of model's parameters to match the real capacitance values.

### 4. `electronicMeasures.m`

This script evaluates the voltage measurements corresponding to different capacitance values from an electronic setup. 

**Features:**
- Bar plots for voltage (mean and standard deviation) for different capacitance values.
- Model fitting.
- Transduction diagram showing the relationship between capacitance and voltage.
- Calibration diagram with uncertainty bands.

## Requirements

- MATLAB R2021b or later.
- RF and Mixed Signal Toolbox for MATLAB.

## Authors

- Carla Finocchiaro
- Emanuele Cannizzo
- Miriam Di Mauro
