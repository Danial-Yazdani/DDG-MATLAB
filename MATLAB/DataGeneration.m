%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
%Author: 
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
% e-mail: X DOT Y AT gmail DOT com
% Copyright notice: (c) 2024 X Y
%**************************************************************************
% function [DDG] = DataGeneration(NewSampleSize,DDG)
DataSample = NaN(NewSampleSize,DDG.NumberOfVariables);
counter=0;
Probability=arrayfun(@(x) x.Weight, DDG.RGC)/sum(arrayfun(@(x) x.Weight, DDG.RGC));
while counter <NewSampleSize
    ChosenID=randsample(DDG.Rng,DDG.RGCNumber, 1, true, Probability);
    Sample=((randn(DDG.Rng,1,DDG.NumberOfVariables).*DDG.RGC(ChosenID).Sigma))*DDG.RGC(ChosenID).RotationMatrix+DDG.RGC(ChosenID).Center;
    if all(Sample >= DDG.MinCoordinate & Sample <= DDG.MaxCoordinate)
       counter=counter+1;
       DataSample(counter,:)=Sample;
    end
end
DDG.Data.Dataset = [DataSample; DDG.Data.Dataset];% Add new samples to the beginning of the dataset
DDG.Data.Dataset = DDG.Data.Dataset(1:end-NewSampleSize, :);% Remove the last NewSampleSize samples to maintain the dataset size