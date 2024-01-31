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
function [DDG] = DataGeneration(NewSampleSize,DDG)
DataSample = NaN(NewSampleSize,DDG.NumberOfVariables);
if DDG.MovingObjects == 1
    DataOriginRandomVector = NaN(NewSampleSize,DDG.NumberOfVariables);
    DataOriginDGCid = NaN(NewSampleSize,1);
end
counter = 0;
Probability = arrayfun(@(x) x.Weight, DDG.DGC)/sum(arrayfun(@(x) x.Weight, DDG.DGC));%Defining the probability of choosing each DGC for generating a data point based on their weight values.
while counter <NewSampleSize
    ChosenID=randsample(DDG.Rng,DDG.DGCNumber, 1, true, Probability);
    RandomVector = randn(DDG.Rng,1,DDG.NumberOfVariables);
    Sample = ((RandomVector.*DDG.DGC(ChosenID).Sigma)*DDG.DGC(ChosenID).RotationMatrix)+DDG.DGC(ChosenID).Center;
    if all(Sample >= DDG.MinCoordinate & Sample <= DDG.MaxCoordinate)
        counter=counter+1;
        DataSample(counter,:) = Sample;
        if DDG.MovingObjects == 1
            DataOriginRandomVector(counter,:) = RandomVector;
            DataOriginDGCid(counter,1) = ChosenID;
        end
    end
end
DDG.Data.Dataset = [DataSample; DDG.Data.Dataset];% Add new samples to the beginning of the dataset
DDG.Data.Dataset = DDG.Data.Dataset(1:DDG.Data.Size, :);% Remove the last NewSampleSize samples to maintain the dataset size
if DDG.MovingObjects == 1
    DDG.Data.Origin.RandomVector  = DataOriginRandomVector;  
    DDG.Data.Origin.DGCid         = DataOriginDGCid;
end