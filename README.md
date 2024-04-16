# Dynamic Dataset Generator (DDG)

### Description
The Dynamic Dataset Generator (DDG) is a MATLAB implementation designed to simulate a wide range of dynamic clustering scenarios. 
This tool facilitates the generation of datasets with heterogeneous, controllable changes, making it ideal for evaluating clustering algorithms in dynamic environments. 
DDG integrates multiple dynamic Gaussian components, which vary in terms of location, scale, and rotation, to reflect realistic and diverse dynamics.

This project is based on the paper:
"Clustering in Dynamic Environments: A Framework for Benchmark Dataset Generation With Heterogeneous Changes" by Danial Yazdani et al., 2024. [Read the paper here](https://arxiv.org/abs/2402.15731v2).

### Prerequisites
- MATLAB (recommended version R2021a or later)

### Setup
1. Clone the repository to your local machine using GitHub or download the ZIP file.

### Usage
To use the DDG, follow these steps:

1. Open the `main.m` file in MATLAB.
2. Run the script to start generating datasets. Configuration options can be adjusted within the script in 'DDGinitialization.m'.

## Features
- **Dynamic Changes**: Simulate datasets with changes in cluster centroids, variances, and number of clusters over time.
- **Customizable Parameters**: Control the degree of dynamism, including the frequency and magnitude of changes.

## License
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

## Authors
- **Danial Yazdani** - *Initial work* - [DanialYazdani](https://github.com/Danial-Yazdani)
