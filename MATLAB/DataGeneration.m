%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
% Author: Danial Yazdani
% Last Edited: November 28, 2025
% Title: Main function of DDG
% --------
% Reference: "Clustering in Dynamic Environments: A Framework for Benchmark
%            Dataset Generation With Heterogeneous Changes"
%            Danial Yazdani et al., 2024. 
%            Available at https://arxiv.org/abs/2402.15731v2
%
% --------
% Description: This function Generates data and updates the dataset.
%
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial.yazdani@gmail.com
% Copyright notice: (c) 2024 Danial Yazdani
%**************************************************************************
function [DDG] = DataGeneration(NewSampleSize,DDG)
% Extract weights into array (faster than arrayfun)
Weights = [DDG.DGC.Weight];
Probability = Weights / sum(Weights); % Normalized probabilities for DGC selection

% Batch select which DGC generates each sample
ChosenIDs = randsample(DDG.Rng, DDG.DGCNumber, NewSampleSize, true, Probability);

% Pre-allocate output
DataSample = zeros(NewSampleSize, DDG.NumberOfVariables);

% Generate all random vectors at once
RandomVectors = randn(DDG.Rng, NewSampleSize, DDG.NumberOfVariables);

% Process each DGC's samples in batch (vectorized)
for dgcId = 1:DDG.DGCNumber
    mask = (ChosenIDs == dgcId);
    nSamples = sum(mask);
    if nSamples > 0
        % Batch generate all samples for this DGC
        DataSample(mask, :) = (RandomVectors(mask, :) .* DDG.DGC(dgcId).Sigma) * DDG.DGC(dgcId).RotationMatrix + DDG.DGC(dgcId).Center;
    end
end

% Update dataset using FIFO approach
DDG.Data.Dataset = [DataSample; DDG.Data.Dataset(1:end-NewSampleSize, :)];
end