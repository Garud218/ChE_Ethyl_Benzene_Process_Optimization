# ChE_Ethyl_Benzene_Process_Optimization

## Project Overview
This repository hosts the detailed implementation and documentation of the **Ethyl Benzene Process Optimization** project, undertaken as a course project for **ChE352** under the supervision of **Prof. Nitin Kaistha**. The project spanned from **March 2025 to April 2025** and focuses on optimizing the ethyl benzene production process through advanced control strategies, optimization techniques, and industrial-scale design validation.

## Objectives
The primary goals of this project are:
- Develop an optimization framework leveraging sensitivity analysis and MATLAB's `fmincon` function to achieve a **47.6% reduction in Total Annual Cost (TAC)**.
- Design and implement dual-control systems using **Ziegler-Nichols** and **Tyreus-Luyben** tuning methods to ensure **99.9% purity** in the output.
- Optimize a **three-column distillation network** with heat integration, resulting in energy savings and a **3-year payback period**.
- Validate the industrial-scale design through comprehensive **material and energy balance calculations**.

## Features
- **Optimization Framework**: Utilizes MATLAB's `fmincon` for nonlinear optimization, incorporating sensitivity analysis to minimize TAC.
- **Control Systems**: Implements dual-control strategies with Ziegler-Nichols and Tyreus-Luyben methods for enhanced process stability and purity.
- **Distillation Network**: Designs a three-column system with heat integration to improve energy efficiency.
- **Validation**: Includes detailed material and energy balance calculations to ensure scalability and practicality for industrial applications.

## Installation and Setup
### Prerequisites
- MATLAB with the **Optimization Toolbox** installed.
- Basic understanding of chemical engineering principles and process optimization.

### Steps to Run
1. **Clone the Repository**:
   ```bash
   git clone https://github.com/Garud218/ChE_Ethyl_Benzene_Process_Optimization.git
   ```
2. **Navigate to the Project Directory**:
   ```bash
   cd ChE_Ethyl_Benzene_Process_Optimization
   ```
3. **Open MATLAB** and load the project workspace.
4. **Execute the Main Scripts**:
   - Run `optimization_script.m` for the TAC reduction analysis.
   - Run `control_design_script.m` for dual-control system simulation.
   - Run `distillation_optimization.m` for the three-column network optimization.
   - Review `balance_calculations.m` for material and energy balance validation.
5. **Analyze Results**: Output files and plots will be generated in the `results` folder.

## Project Structure
- `src/`: Contains all MATLAB scripts and functions.
  - `optimization_script.m`: Main script for TAC optimization.
  - `control_design_script.m`: Implements Ziegler-Nichols and Tyreus-Luyben control systems.
  - `distillation_optimization.m`: Optimizes the three-column distillation network.
  - `balance_calculations.m`: Performs material and energy balance calculations.
- `data/`: Input data files and parameters.
- `results/`: Output files, including optimization results and validation data.
- `docs/`: Additional documentation and reports.
- `README.md`: This file.

## Results
- **TAC Reduction**: Successfully reduced TAC by **47.6%** through the optimization framework.
- **Purity Achievement**: Attained **99.9% purity** using the dual-control systems.
- **Energy Efficiency**: Heat integration in the distillation network led to a **3-year payback period** with significant energy savings.
- **Validation**: Industrial-scale design validated with accurate material and energy balance calculations.

## Contributing
We welcome contributions to enhance this project! To contribute:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -m "Add new feature"`).
4. Push to the branch (`git push origin feature-branch`).
5. Open a pull request with a detailed description of your changes.

## Acknowledgments
- **Prof. Nitin Kaistha** for guidance and mentorship.
- The ChE352 course team for providing the opportunity to work on this project.
- MATLAB for providing the necessary tools for optimization and simulation.

## Contact
For any questions or suggestions, please open an issue on the GitHub repository or contact the project maintainer at [Garud218](https://github.com/Garud218).
