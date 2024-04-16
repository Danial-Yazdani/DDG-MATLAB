%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
% Author: Danial Yazdani
% Last Edited: January 31, 2024
% Title: Main function of DDG
% --------
% Reference: "Clustering in Dynamic Environments: A Framework for Benchmark
%            Dataset Generation With Heterogeneous Changes"
%            Danial Yazdani et al., 2024. 
%            Available at https://arxiv.org/abs/2402.15731v2
%
% --------
% Description:
% This MATLAB script implements the Dynamic Dataset Generator (DDG), which
% simulates a wide range of dynamic clustering scenarios. DDG allows for the
% generation of datasets with heterogeneous, controllable changes, making it
% ideal for evaluating clustering algorithms in dynamic environments. The generator
% integrates multiple dynamic Gaussian components, which vary in terms of location,
% scale, and rotation, reflecting realistic and diverse dynamics.
%
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial.yazdani@gmail.com
% Copyright notice: (c) 2024 Danial Yazdani
%**************************************************************************
close all;clear all;clc; %#ok<CLALL>
%% Preparation and initialization
clear DDG;
DDG = DDGinitialization;
rng('shuffle');%Set a random seed for the optimizer
%% Start optimization/clustering method code from here


%Your optimization/clustering method code%

%OfflinePerformance = mean(DDG.CurrentBestSolutionValue);% The performance of an algorithm in a single run.