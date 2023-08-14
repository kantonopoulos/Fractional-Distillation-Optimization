# Fractional Distillation Optimization
***
## Table of Contents
1. [Introduction](#introduction)
2. [What can I do with the Rheological Analysis Toolbox?](#what-can-i-do-with-the-rheological-analysis-toolbox?)
3. [License](#license)
4. [Installation](#installation)

## Overview
The Fractional Distillation Optimization is a powerful in-house algorithm designed to determine the optimal design characteristics of a fractional distillation column. The algorithm focuses on refining key parameters such as reflux ratio, number of trays, feed location, and more. This code has been developed using the Lewis-Matheson method, material and energy balances, and vapor-liquid equilibrium (VLE) data. It has been specifically tailored for fractional distillation scenarios, with a primary application in the separation of an ethanol-water mixture.

## Features
* **Optimization of Design Parameters**: The algorithm optimizes essential design parameters of a distillation column, including the reflux ratio, number of trays, and feed location. This enables the creation of an efficient column that meets specific separation requirements.
* **Comprehensive Design Calculations**: The code implements material and energy balances, as well as vapor-liquid equilibrium (VLE) data, to ensure accurate and comprehensive design calculations. This results in a distillation column that operates reliably and effectively.
* **Physical Process Modeling**: Drawing inspiration from the work of Assael and Maggiliotou in "Physical Processes: A Calculation Introduction," the algorithm embodies a solid understanding of the underlying physical processes involved in fractional distillation.
* **Economical Assessment**: The algorithm evaluates the economic viability of the designed distillation column. It considers costs associated with utilities (steam and cooling water), annual revenue loss, and overall operating costs over the expected lifetime of the column.
* **Scenario-Based Approach**: The algorithm allows for the exploration of multiple scenarios by varying the mass fractions in the feed, distillate, and bottoms. This flexibility enables users to identify the best design under different conditions.

## Contribution
We acknowledge the contribution of Assael and Maggiliotou for their valuable insights and foundational work in the field of physical processes, which have greatly influenced the development of this algorithm.

## How to Use
1. Make sure you have Python installed on your system.
2. Install the required libraries (NumPy and Pandas) using the command:
```
pip install numpy pandas
```
3. Copy and paste the provided code into a **`.py`** file (e.g., **`fractional_distillation_optimization.py`**).
Run the script using a Python interpreter:
```
python fractional_distillation_optimization.py
```
The algorithm will generate an Excel file named **`distcol.xlsx`** containing optimal design parameters for various scenarios.

## Important Notes
* This algorithm is specifically designed for educational purposes and optimized for ethanol-water mixture separation. It might require adjustments for different compound mixtures or specific constraints.
* We recommend adding error handling and explanatory comments for production use, enhancing the robustness and clarity of the code.
* The algorithm relies on accurate physical properties and models. Ensure the relevancy and accuracy of these properties for your application.
* The optimization process might take some time depending on the scenarios and mass fractions considered. Modify the mass fraction ranges and step sizes to achieve the desired balance between accuracy and computation time.
