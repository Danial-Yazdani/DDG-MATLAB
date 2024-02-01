%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
% Author: X Y
%Last Edited: January 31, 2024
%Title: Main function of DDG
% --------
%Refrence: "Clustering in Dynamic Environments: A Framework for Benchmark
%          Dataset Generation With Heterogeneous Changes"
%
% --------
% Description: This function Generates data and updates the dataset.
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: X Y
% e-mail: X DOT Y AT something DOT com
% Copyright notice: (c) 2024 X Y
%**************************************************************************
function [DDG] = MovingObjects(DDG,DGCid)
for ii=1 : DDG.Data.Size
    if DDG.Data.Origin.DGCid(ii,1) == DGCid
        if rand(DDG.Rng) > DDG.MovingObjectsResamplingLikelihood
            DDG.Data.Dataset(ii,:) = ((DDG.Data.Origin.RandomVector(ii,:).*DDG.DGC(DGCid).Sigma)*DDG.DGC(DGCid).RotationMatrix)+DDG.DGC(DGCid).Center;
        else
            DDG.Data.Origin.RandomVector(ii,:) = randn(DDG.Rng,1,DDG.NumberOfVariables);
            DDG.Data.Dataset(ii,:) = ((DDG.Data.Origin.RandomVector(ii,:).*DDG.DGC(DGCid).Sigma)*DDG.DGC(DGCid).RotationMatrix)+DDG.DGC(DGCid).Center;
        end
    end
end