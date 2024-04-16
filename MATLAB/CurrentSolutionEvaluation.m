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
% Description: This function evaluates a single clustering solution, employing
% the same objective function and representation as ClusteringEvaluation(.).
% Its primary use is to re-evaluate the current best solution following
% updates to the dataset, specifically for performance measurement purposes.
% This separate implementation is designed to avoid recursion, which would 
% occur if ClusteringEvaluation(.) called itself. Importantly, as this 
% function is used solely for performance measurement, its usage is not 
% counted in the total number of function evaluations.
%
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial.yazdani@gmail.com
% Copyright notice: (c) 2024 Danial Yazdani
%**************************************************************************
function [result] = CurrentSolutionEvaluation(x,DDG)
ClusterCenterPosition = reshape(x', [DDG.NumberOfVariables, DDG.ClusterNumber])';
Distances = pdist2(DDG.Data.Dataset, ClusterCenterPosition,'euclidean');
[~, closestClusterIndices] = min(Distances, [], 2);
selectedDistances = diag(Distances(:, closestClusterIndices));
result = sum(selectedDistances);% Sum of intra-cluster distances
end