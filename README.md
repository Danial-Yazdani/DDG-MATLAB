# Dynamic Dataset Generator (DDG) - MATLAB Implementation

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a+-orange.svg)](https://www.mathworks.com/products/matlab.html)

## Overview

The **Dynamic Dataset Generator (DDG)** is a comprehensive MATLAB framework for generating dynamic datasets with controllable characteristics. It is specifically designed for benchmarking and evaluating clustering algorithms in dynamic environments where data distributions evolve over time.

DDG simulates realistic dynamic scenarios using **Dynamic Gaussian Components (DGCs)**, which can vary in location, scale, rotation, and weight, enabling the creation of diverse and challenging benchmark datasets.

**Paper**: [Clustering in Dynamic Environments: A Framework for Benchmark Dataset Generation With Heterogeneous Changes](https://arxiv.org/abs/2402.15731v2)  
**Published in**: GECCO 2024 (Genetic and Evolutionary Computation Conference)  
*Danial Yazdani et al., 2024*

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [File Structure](#file-structure)
- [Configuration Guide](#configuration-guide)
- [API Reference](#api-reference)
- [Dynamic Change Types](#dynamic-change-types)
- [Examples](#examples)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Features

- **Multiple Dynamic Gaussian Components (DGCs)**: Generate multimodal data distributions
- **Gradual Local Changes**: Smooth, correlated changes to individual DGC parameters
- **Severe Global Changes**: Abrupt changes affecting all DGCs simultaneously
- **Configurable Dynamics**: Control over number of DGCs, variables, and clusters
- **Rotation Support**: Full rotation matrix control for each DGC
- **Boundary Control**: Reflect method ensures parameters stay within valid ranges
- **Performance Tracking**: Built-in evaluation and performance measurement
- **Reproducible Results**: Seed-based random number generation

---

## Installation

### Prerequisites
- MATLAB R2021a or later (earlier versions may work but are not tested)
- Statistics and Machine Learning Toolbox (for `betarnd` function)

### Setup
1. Clone or download this repository:
   ```bash
   git clone https://github.com/Danial-Yazdani/DDG-MATLAB.git
   ```
2. Add the `MATLAB` folder to your MATLAB path:
   ```matlab
   addpath('path/to/DDG-MATLAB-main/MATLAB')
   ```

---

## Quick Start

```matlab
% Initialize DDG with default settings
DDG = DDGinitialization();

% Your clustering algorithm goes here
% Example: Evaluate a random solution
solution = rand(1, DDG.ClusterNumber * DDG.NumberOfVariables) * 140 - 70;
[fitness, DDG] = ClusteringEvaluation(solution, DDG);

% Access the generated dataset
dataset = DDG.Data.Dataset;  % Size: [DDG.Data.Size x DDG.NumberOfVariables]
```

---

## File Structure

| File | Description |
|------|-------------|
| `main.m` | Entry point - initialize DDG and run your algorithm |
| `DDGinitialization.m` | Configure all DDG parameters and initialize DGCs |
| `DataGeneration.m` | Generate data samples from DGCs |
| `ClusteringEvaluation.m` | Evaluate clustering solutions and trigger changes |
| `EnvironmentalChangeGenerator.m` | Apply environmental changes to DGCs |
| `CurrentSolutionEvaluation.m` | Re-evaluate solutions after dataset updates |

---

## Configuration Guide

All parameters are configured in `DDGinitialization.m`. Key parameters include:

### Basic Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DDG.Seed` | 2151 | Random seed for reproducibility |
| `DDG.MaxEvals` | 500000 | Maximum function evaluations |
| `DDG.NumberOfVariables` | 2 | Dimensionality of data |
| `DDG.DGCNumber` | 7 | Number of Dynamic Gaussian Components |
| `DDG.ClusterNumber` | 5 | Number of clusters |

### Bounds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DDG.MinCoordinate` / `MaxCoordinate` | -70 / 70 | DGC center bounds |
| `DDG.MinSigma` / `MaxSigma` | 7 / 20 | Standard deviation bounds |
| `DDG.MinWeight` / `MaxWeight` | 1 / 3 | DGC weight bounds |
| `DDG.MinAngle` / `MaxAngle` | -π / π | Rotation angle bounds |

### Change Severity (Local/Gradual)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DDG.LocalShiftSeverityRange` | [0.1, 0.2] | Center shift magnitude |
| `DDG.RelocationCorrelationRange` | [0.99, 0.995] | Movement correlation (ρ) |
| `DDG.LocalSigmaSeverityRange` | [0.05, 0.1] | Sigma change magnitude |
| `DDG.LocalWeightSeverityRange` | [0.02, 0.05] | Weight change magnitude |
| `DDG.LocalRotationSeverityRange` | [π/360, π/180] | Rotation change magnitude |

### Change Severity (Global/Severe)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DDG.GlobalShiftSeverityValue` | 10 | Global center shift |
| `DDG.GlobalSigmaSeverityValue` | 5 | Global sigma change |
| `DDG.GlobalWeightSeverityValue` | 0.5 | Global weight change |
| `DDG.GlobalAngleSeverityValue` | π/4 | Global rotation change |
| `DDG.GlobalSeverityControl` | 0.1 | Beta distribution parameter (α=β) |

### Change Likelihoods

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DDG.LocalTemporalSeverityRange` | [0.05, 0.1] | Probability of local change per DGC |
| `DDG.GlobalChangeLikelihood` | 0.0001 | Probability of global change |
| `DDG.DGCNumberChangeLikelihood` | 0.0001 | Probability of DGC count change |
| `DDG.VariableNumberChangeLikelihood` | 0.0001 | Probability of dimension change |
| `DDG.ClusterNumberChangeLikelihood` | 0.0001 | Probability of cluster count change |

### Dataset Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DDG.Data.Size` | 1000 | Dataset size |
| `DDG.Data.FrequentSamplingLikelihood` | 0.1 | Incremental sampling probability |
| `DDG.Data.IncrementalSamplingSize` | 50 (5%) | Samples added per incremental update |

---

## API Reference

### `DDGinitialization()`
Initializes and returns the DDG structure with all parameters and initial dataset.

**Returns**: `DDG` structure

---

### `ClusteringEvaluation(X, DDG)`
Evaluates clustering solutions and triggers environmental changes.

**Parameters**:
- `X`: Solution matrix [n_solutions × (n_clusters × n_variables)]
- `DDG`: DDG structure

**Returns**:
- `result`: Fitness values (sum of intra-cluster distances)
- `DDG`: Updated DDG structure

---

### `DataGeneration(NewSampleSize, DDG)`
Generates new data samples from current DGCs.

**Parameters**:
- `NewSampleSize`: Number of samples to generate
- `DDG`: DDG structure

**Returns**: Updated `DDG` with new samples added to dataset

---

### `EnvironmentalChangeGenerator(DDG, ChangeCode)`
Applies environmental changes based on change code.

**Change Codes**:
- `> 0`: Local change for DGC with ID = ChangeCode
- `0`: Global severe change for all DGCs
- `-1`: Change number of DGCs
- `-2`: Change number of variables
- `-3`: Change number of clusters

---

## Dynamic Change Types

### 1. Gradual Local Changes (per DGC)
- **Center Relocation**: Correlated random walk (Equations 5-6)
- **Sigma Changes**: Direction-based incremental changes (Equation 7)
- **Weight Changes**: Direction-based incremental changes (Equation 8)
- **Rotation Changes**: Angle adjustments (Equation 9)

### 2. Severe Global Changes (all DGCs)
- Heavy-tail Beta distribution for significant shifts (Equations 10-13)
- Applied with probability `DDG.GlobalChangeLikelihood`

### 3. Structural Changes
- Number of DGCs can increase/decrease (Equation 14)
- Number of variables can change (Equation 15)
- Number of clusters can change (Equation 16)

---

## Examples

### Example 1: Basic Usage
```matlab
% Initialize
DDG = DDGinitialization();

% Run optimization loop
for iter = 1:1000
    % Generate candidate solutions
    solutions = rand(10, DDG.ClusterNumber * DDG.NumberOfVariables) * 140 - 70;
    
    % Evaluate
    [fitness, DDG] = ClusteringEvaluation(solutions, DDG);
    
    % Check termination
    if DDG.FE >= DDG.MaxEvals
        break;
    end
end

% Results
fprintf('Best fitness: %.4f\n', DDG.CurrentBestSolutionValue);
```

### Example 2: Static Environment (No Changes)
```matlab
DDG = DDGinitialization();

% Disable all changes
DDG.GlobalChangeLikelihood = 0;
DDG.DGCNumberChangeLikelihood = 0;
DDG.VariableNumberChangeLikelihood = 0;
DDG.ClusterNumberChangeLikelihood = 0;
DDG.Data.FrequentSamplingLikelihood = 0;
for i = 1:DDG.DGCNumber
    DDG.DGC(i).LocalChangeLikelihood = 0;
end
```

### Example 3: High-Frequency Changes
```matlab
DDG = DDGinitialization();

% Increase change frequency
DDG.GlobalChangeLikelihood = 0.001;
DDG.LocalTemporalSeverityRange = [0.2, 0.3];
```

---

## Citation

If you use DDG in your research, please cite:

```bibtex
@inproceedings{yazdani2024clustering,
  title={Clustering in Dynamic Environments: A Framework for Benchmark Dataset Generation With Heterogeneous Changes},
  author={Yazdani, Danial and Branke, J{\"u}rgen and Meghdadi, Amir Hossein and Omidvar, Mohammad Nabi and Stoean, Catalin and Stoean, Ruxandra and Gandomi, Amir H and Shi, Yuhui},
  booktitle={Proceedings of the Genetic and Evolutionary Computation Conference (GECCO)},
  year={2024},
  organization={ACM}
}
```

**arXiv preprint**: [arXiv:2402.15731](https://arxiv.org/abs/2402.15731v2)

---

## License

This project is licensed under the **GNU General Public License v3.0** - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Danial Yazdani**  
Email: danial.yazdani@gmail.com  
GitHub: [Danial-Yazdani](https://github.com/Danial-Yazdani)

---

## Acknowledgments

- Thanks to all contributors and users who provide feedback
- Python implementation also available: [DDG-Python](https://github.com/Danial-Yazdani/DDG-Python)
