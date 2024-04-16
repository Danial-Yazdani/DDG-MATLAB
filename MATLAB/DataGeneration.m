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
DataSample = NaN(NewSampleSize,DDG.NumberOfVariables);
counter = 0;
Probability = arrayfun(@(x) x.Weight, DDG.DGC)/sum(arrayfun(@(x) x.Weight, DDG.DGC));%Defining the probability of choosing each DGC for generating a data point based on their weight values.
while counter <NewSampleSize
    ChosenID=randsample(DDG.Rng,DDG.DGCNumber, 1, true, Probability);
    RandomVector = randn(DDG.Rng,1,DDG.NumberOfVariables);
    Sample = ((RandomVector .* DDG.DGC(ChosenID).Sigma) * DDG.DGC(ChosenID).RotationMatrix) + DDG.DGC(ChosenID).Center;
    %     if all(Sample >= DDG.MinCoordinate & Sample <= DDG.MaxCoordinate) % Activate this IF if you want to keep the data points inside the boundaries of the MEAN positions
    counter=counter+1;
    DataSample(counter,:) = Sample;
    %     end
end
DDG.Data.Dataset = [DataSample; DDG.Data.Dataset];% Add new samples to the beginning of the dataset
DDG.Data.Dataset = DDG.Data.Dataset(1:DDG.Data.Size, :);% Remove the last NewSampleSize samples to maintain the dataset size
end