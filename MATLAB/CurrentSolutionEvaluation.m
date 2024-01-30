%***********************************DDG ver 2.00**************************
%Author: Danial Yazdani
%Last Edited: December 14, 2023
%Title:
% --------
%Refrence:
%
%
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial DOT yazdani AT gmail DOT com
% Copyright notice: (c) 2023 Danial Yazdani
%**************************************************************************
function [result] = CurrentSolutionEvaluation(x,DDG)
ClusterCenterPosition = reshape(x', [DDG.NumberOfVariables, DDG.ClusterNumber])';
Distances = pdist2(DDG.Data.Dataset, ClusterCenterPosition,'euclidean');
[~, closestClusterIndices] = min(Distances, [], 2);
selectedDistances = diag(Distances(:, closestClusterIndices));
result = sum(selectedDistances);% Sum of intra-cluster distances
end