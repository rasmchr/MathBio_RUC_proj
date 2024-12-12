# MathBio_RUC_proj
# Penicillin Production Modeling Project

## Overview
This repository contains a collection of Python, MATLAB, and Maple scripts designed to perform computational analysis and visualizations for the mathematical modeling of Penicillin G production using industrial-scale fermentation data. The scripts facilitate model development, parameter estimation, numerical simulations, and visual representation of results, enabling better understanding and optimization of the fermentation process.

## Project Introduction
Since Alexander Flemming’s discovery of Penicillin’s antibiotic properties in 1929, significant progress has been made in improving its yield and production methods. Over the past 70 years, Penicillium chrysogenum (P. chrysogenum) has been used extensively in large-scale production, with mathematical modeling playing a critical role in optimizing the process. Mathematical models describe relationships between key variables, such as biomass, substrate, and product, often through ordinary differential equations (ODEs).

This project aims to construct a mathematical model that captures the dynamics observed in industrial-scale fermentation processes of Penicillin G. By leveraging industrial-scale data and revisiting foundational models, this work highlights the relevance of simple, low-dimensional models in providing reliable predictions. These models are particularly useful during the initial stages of fermentation process development, reducing costs and enabling early-stage decision-making.

The research question guiding this project is:
> How is a low-dimensional, coarse-grain model able to reflect the dynamics of a real industrial-scale fermentation process of Penicillium chrysogenum?

## Repository Contents
This repository is structured as follows:

### 1. **Python Scripts**
- Scripts for numerical simulations of ODE-based models.
- Data analysis and visualization utilities.
- Parameter estimation using optimization techniques.

### 2. **MATLAB Scripts**
- Implementation of models in MATLAB for verification and cross-platform comparison.
- Tools for sensitivity analysis and control strategy design.

### 3. **Maple Scripts**
- Symbolic computation for deriving model equations and simplifying expressions.
- Analytical exploration of system dynamics.

## Key Features
- **Cross-Platform Compatibility**: Use of Python, MATLAB, and Maple to leverage the strengths of each platform.
- **Data-Driven Modeling**: Incorporates industrial-scale data for realistic simulation and validation.
- **Comprehensive Analysis**: Supports numerical, analytical, and visualization tasks to facilitate understanding of the fermentation process.
- **Revisiting Foundational Models**: Builds upon and modernizes existing models for contemporary use cases.

## Getting Started

### Prerequisites
- Install Python (3.7 or higher) with the required libraries (see `requirements.txt`).
- MATLAB R2021a or higher with the necessary toolboxes.
- Maple 2020 or higher.

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/penicillin-modeling.git
   ```
2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Ensure MATLAB and Maple are installed and properly configured.

### Usage
- Use Python scripts for primary data analysis and visualization.
- Use MATLAB scripts for sensitivity analysis and advanced control design.
- Use Maple scripts for symbolic computations and analytical explorations.

## References
1. Bajpai & Reuß (1980): Foundational model of Penicillin G production.
2. Industrial-scale data source: Reference [7].
3. Additional references: [10], [15], [25], and [13] from the project’s introduction.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments
Special thanks to the authors of the referenced works and contributors to the modeling community for their foundational contributions. Their efforts have inspired this project’s approach and methodology.

---
For any questions or contributions, feel free to contact the repository maintainer or submit an issue through the GitHub repository.

