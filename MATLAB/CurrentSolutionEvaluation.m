%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
%Author: 
%Last Edited: January 31, 2024
%Title: Main function of DDG
% --------
%Refrence: "Clustering in Dynamic Environments: A Framework for Benchmark 
%          Dataset Generation With Heterogeneous Changes"
%
% --------
% Description: This function evaluates a single clustering solution, utilizing
% the same objective function and representation as ClusteringEvaluation(.).
% It is primarily used for re-evaluating the current best solution after
% updates to the dataset. This implementation serves to prevent recursion 
% that would occur if ClusteringEvaluation(.) were to call itself. 
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: X Y
% e-mail: X DOT Y AT gmail DOT com
% Copyright notice: (c) 2024 X Y
%**************************************************************************
function [result] = CurrentSolutionEvaluation(x,DDG)
ClusterCenterPosition = reshape(x', [DDG.NumberOfVariables, DDG.ClusterNumber])';
Distances = pdist2(DDG.Data.Dataset, ClusterCenterPosition,'euclidean');
[~, closestClusterIndices] = min(Distances, [], 2);
selectedDistances = diag(Distances(:, closestClusterIndices));
result = sum(selectedDistances);% Sum of intra-cluster distances
end