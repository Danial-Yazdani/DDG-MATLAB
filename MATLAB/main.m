%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
% Author: X Y
%Last Edited: January 31, 2024
%Title: Main function of DDG
% --------
%Refrence: "Clustering in Dynamic Environments: A Framework for Benchmark
%          Dataset Generation With Heterogeneous Changes"
%
%
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: X Y
% e-mail: X DOT Y AT something DOT com
% Copyright notice: (c) 2024 X Y
%**************************************************************************
close all;clear all;clc; %#ok<CLALL>
%% Preparation and initialization
clear DDG;
DDG = DDGinitialization;
rng('shuffle');%Set a random seed for the optimizer
%% Start optimization/clustering method code from here


%Your optimization/clustering method code%

%OfflinePerformance = mean(DDG.CurrentBestSolutionValue);% The performance of an algorithm in a single run.