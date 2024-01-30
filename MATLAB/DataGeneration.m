function [DDG] = DataGeneration(NewSampleSize,DDG)
DataSample = NaN(NewSampleSize,DDG.NumberOfVariables);
counter=0;
weights=arrayfun(@(x) x.Weight, DDG.RGC)/sum(arrayfun(@(x) x.Weight, DDG.RGC));
while counter <NewSampleSize
    ChosenID=randsample(DDG.Rng,DDG.RGCNumber, 1, true, weights);
    Sample=((randn(DDG.Rng,1,DDG.NumberOfVariables).*DDG.RGC(ChosenID).Sigma))*DDG.RGC(ChosenID).RotationMatrix+DDG.RGC(ChosenID).Center;
    if all(Sample >= DDG.MinCoordinate & Sample <= DDG.MaxCoordinate)
       counter=counter+1;
       DataSample(counter,:)=Sample;
    end
end
DDG.Data.Dataset = [DataSample; DDG.Data.Dataset];% Add new samples to the beginning of the dataset
DDG.Data.Dataset = DDG.Data.Dataset(1:end-NewSampleSize, :);% Remove the last NewSampleSize samples to maintain the dataset size