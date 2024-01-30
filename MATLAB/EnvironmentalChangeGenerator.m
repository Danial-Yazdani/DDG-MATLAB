%***********************************DDG ver 1.00**************************
%Author: Danial Yazdani
%Last Edited: December 14, 2023
%Title: Environmental change in the dynamic landscape generator
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
function [DDG] = EnvironmentalChangeGenerator(DDG,ChangeCode)
%The default dynamic used for generating environmental changes is random.
if ChangeCode>0
    %% Gradual local change for the DGC whose ID is equal to ChangeCode
    RGCid = ChangeCode;
    % RGC center update
    tmp = randn(DDG.Rng,1,DDG.NumberOfVariables);
    RandomDirection = tmp/sqrt(sum(tmp.^2));
    SummedVector=(((1-DDG.RGC(RGCid).ShiftCorrelationFactor)*RandomDirection)+ (DDG.RGC(RGCid).ShiftCorrelationFactor*DDG.RGC(RGCid).PreviousShiftDirection));
    RelocationDirection = SummedVector/ sqrt(sum(SummedVector.^2));% Normalized to create a unit vector
    UpdatedRGCsPosition  = DDG.RGC(RGCid).Center + (RelocationDirection*(abs(randn(DDG.Rng))*DDG.RGC(RGCid).ShiftSeverity));
    % Boundary check for DGC center
    tmp = UpdatedRGCsPosition > DDG.MaxCoordinate;
    UpdatedRGCsPosition(tmp) = (2*DDG.MaxCoordinate)- UpdatedRGCsPosition(tmp);
    tmp = UpdatedRGCsPosition < DDG.MinCoordinate;
    UpdatedRGCsPosition(tmp) = (2*DDG.MinCoordinate)- UpdatedRGCsPosition(tmp);
    RelocationVectore = UpdatedRGCsPosition-DDG.RGC(RGCid).Center;
    DDG.RGC(RGCid).PreviousShiftDirection = (RelocationVectore)./ sqrt(sum(RelocationVectore.^2));%Storing the last relocation direction
    DDG.RGC(RGCid).Center = UpdatedRGCsPosition;
    % Weight update
    if rand(DDG.Rng)<DDG.RGC(RGCid).DirectionChangeProbabolity
        DDG.RGC(RGCid).WeightDirection = -DDG.RGC(RGCid).WeightDirection;%Update heigh change direction based on the probability
    end
    UpdatedWeight = DDG.RGC(RGCid).Weight + (abs(randn(DDG.Rng))*DDG.RGC(RGCid).WeightSeverity*DDG.RGC(RGCid).WeightDirection);
    % Boundary check for DGC height
    tmp = UpdatedWeight > DDG.MaxWeight;
    DDG.RGC(RGCid).WeightDirection(tmp) = -DDG.RGC(RGCid).WeightDirection(tmp); % Changing direction if reached the boundary
    UpdatedWeight(tmp) = (2*DDG.MaxWeight)- UpdatedWeight(tmp);
    tmp = UpdatedWeight < DDG.MinWeight;
    DDG.RGC(RGCid).WeightDirection(tmp) = -DDG.RGC(RGCid).WeightDirection(tmp); % Changing direction if reached the boundary
    UpdatedWeight(tmp) = (2*DDG.MinWeight)- UpdatedWeight(tmp);
    DDG.RGC(RGCid).Weight = UpdatedWeight;
    % Sigma update
    if DDG.Conditioning == 0%Keeps condition number at 1
        if rand(DDG.Rng)<DDG.RGC(RGCid).DirectionChangeProbabolity
            DDG.RGC(RGCid).SigmaDirection = -DDG.RGC(RGCid).SigmaDirection;
        end
        UpdatedSigma = DDG.RGC(RGCid).Sigma + ( (ones(1,DDG.NumberOfVariables).*abs(randn(DDG.Rng))) .* DDG.RGC(RGCid).SigmaSeverity .* DDG.RGC(RGCid).SigmaDirection);
    else
        invertFlags = rand(DDG.Rng,1,DDG.NumberOfVariables) < DDG.RGC(RGCid).DirectionChangeProbabolity;%Update the elements of sigma change direction based on the probability
        DDG.RGC(RGCid).SigmaDirection(invertFlags) = -DDG.RGC(RGCid).SigmaDirection(invertFlags);
        UpdatedSigma = DDG.RGC(RGCid).Sigma + (abs(randn(DDG.Rng,1,DDG.NumberOfVariables))*DDG.RGC(RGCid).SigmaSeverity.*DDG.RGC(RGCid).SigmaDirection);
    end
    % Boundary check for DGC sigma
    tmp = UpdatedSigma > DDG.MaxSigma;
    DDG.RGC(RGCid).SigmaDirection(tmp) = -DDG.RGC(RGCid).SigmaDirection(tmp); % Changing direction if reached the boundary
    UpdatedSigma(tmp) = (2*DDG.MaxSigma)- UpdatedSigma(tmp);
    tmp = UpdatedSigma < DDG.MinSigma;
    DDG.RGC(RGCid).SigmaDirection(tmp) = -DDG.RGC(RGCid).SigmaDirection(tmp); % Changing direction if reached the boundary
    UpdatedSigma(tmp) = (2*DDG.MinSigma)- UpdatedSigma(tmp);
    DDG.RGC(RGCid).Sigma = UpdatedSigma;
    % Angle update
    if DDG.Rotation==1
        invertFlags = triu((rand(DDG.Rng,DDG.NumberOfVariables,DDG.NumberOfVariables) < DDG.RGC(RGCid).DirectionChangeProbabolity),1);
        DDG.RGC(RGCid).RotationDirection(invertFlags) = -DDG.RGC(RGCid).RotationDirection(invertFlags);
        UpdatedAngles = DDG.RGC(RGCid).ThetaMatrix + (triu(abs(randn(DDG.Rng,DDG.NumberOfVariables,DDG.NumberOfVariables)),1).*DDG.RGC(RGCid).RotationSeverity.*DDG.RGC(RGCid).RotationDirection);
        % Boundary check for angles
        tmp = UpdatedAngles > DDG.MaxAngle;
        DDG.RGC(RGCid).RotationDirection(tmp) = -DDG.RGC(RGCid).RotationDirection(tmp); % Changing direction if reached the boundary
        UpdatedAngles(tmp) = (2*DDG.MaxAngle)- UpdatedAngles(tmp);
        tmp = UpdatedAngles < DDG.MinAngle;
        DDG.RGC(RGCid).RotationDirection(tmp) = -DDG.RGC(RGCid).RotationDirection(tmp); % Changing direction if reached the boundary
        UpdatedAngles(tmp) = (2*DDG.MinAngle)- UpdatedAngles(tmp);
        DDG.RGC(RGCid).ThetaMatrix = UpdatedAngles;
        [DDG.RGC(RGCid).RotationMatrix] = Rotation(DDG.RGC(RGCid).ThetaMatrix,DDG.NumberOfVariables);
    end
elseif ChangeCode==0 % Severe global change for all DGCs
    RandStream.setGlobalStream(DDG.Rng); % Set the global RandStream to DDG.Rng since betarand does not accpet DDG.Rng as RandStream.
    for RGCid=1 : DDG.RGCNumber
        % RGC center update
        tmp = randn(DDG.Rng,1,DDG.NumberOfVariables);
        RandomDirection = tmp/sqrt(sum(tmp.^2));
        UpdatedRGCsPosition  = DDG.RGC(RGCid).Center + (RandomDirection*DDG.GlobalShiftSeverityValue*(2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl)-1));
        % Boundary check for DGC center
        tmp = UpdatedRGCsPosition > DDG.MaxCoordinate;
        UpdatedRGCsPosition(tmp) = (2*DDG.MaxCoordinate)- UpdatedRGCsPosition(tmp);
        tmp = UpdatedRGCsPosition < DDG.MinCoordinate;
        UpdatedRGCsPosition(tmp) = (2*DDG.MinCoordinate)- UpdatedRGCsPosition(tmp);
        DDG.RGC(RGCid).Center = UpdatedRGCsPosition;
        % Weight update
        UpdatedWeight = DDG.RGC(RGCid).Weight + (DDG.GlobalWeightSeverityValue*(2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl)-1));
        % Boundary check for DGC height
        tmp = UpdatedWeight > DDG.MaxWeight;
        UpdatedWeight(tmp) = (2*DDG.MaxWeight)- UpdatedWeight(tmp);
        tmp = UpdatedWeight < DDG.MinWeight;
        UpdatedWeight(tmp) = (2*DDG.MinWeight)- UpdatedWeight(tmp);
        DDG.RGC(RGCid).Weight = UpdatedWeight;
        % Sigma update
        if DDG.Conditioning == 0%Keeps condition number at 1
            UpdatedSigma = DDG.RGC(RGCid).Sigma + ( (ones(1,DDG.NumberOfVariables).* (2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl)-1)) .* DDG.GlobalSigmaSeverityValue);
        else
            UpdatedSigma = DDG.RGC(RGCid).Sigma + ( (2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl,1,DDG.NumberOfVariables)-1) .* DDG.GlobalSigmaSeverityValue);
        end
        % Boundary check for DGC sigma
        tmp = UpdatedSigma > DDG.MaxSigma;
        UpdatedSigma(tmp) = (2*DDG.MaxSigma)- UpdatedSigma(tmp);
        tmp = UpdatedSigma < DDG.MinSigma;
        UpdatedSigma(tmp) = (2*DDG.MinSigma)- UpdatedSigma(tmp);
        DDG.RGC(RGCid).Sigma = UpdatedSigma;
        % Angle update
        if DDG.Rotation==1
            UpdatedAngles = DDG.RGC(RGCid).ThetaMatrix + (DDG.GlobalAngleSeverityValue*triu((2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl,DDG.NumberOfVariables,DDG.NumberOfVariables)-1),1));
            % Boundary check for angles
            tmp = UpdatedAngles > DDG.MaxAngle;
            UpdatedAngles(tmp) = (2*DDG.MaxAngle)- UpdatedAngles(tmp);
            tmp = UpdatedAngles < DDG.MinAngle;
            UpdatedAngles(tmp) = (2*DDG.MinAngle)- UpdatedAngles(tmp);
            DDG.RGC(RGCid).ThetaMatrix = UpdatedAngles;
            [DDG.RGC(RGCid).RotationMatrix] = Rotation(DDG.RGC(RGCid).ThetaMatrix,DDG.NumberOfVariables);
        end
    end
    rng('shuffle');%Shuffle the global randstream
elseif ChangeCode==-1 % Change in the number of DGCs
    UpdatedRGCNumber = DDG.RGCNumber + sign(randn(DDG.Rng))*DDG.RGCNumberChangeSeverity;
    % Boundary check for DGC sigma
    tmp = UpdatedRGCNumber > DDG.MaxRGCNumber;
    UpdatedRGCNumber(tmp) = (2*DDG.MaxRGCNumber)- UpdatedRGCNumber(tmp);
    tmp = UpdatedRGCNumber < DDG.MinRGCNumber;
    UpdatedRGCNumber(tmp) = (2*DDG.MinRGCNumber)- UpdatedRGCNumber(tmp);
    if UpdatedRGCNumber<DDG.RGCNumber
        RGCsToBeRemoved = randperm(DDG.Rng,DDG.RGCNumber,abs(DDG.RGCNumber-UpdatedRGCNumber));
        DDG.RGC(RGCsToBeRemoved)=[];% Removing randomly chosen DGCs
    elseif UpdatedRGCNumber>DDG.RGCNumber
        for ii=DDG.RGCNumber+1 : UpdatedRGCNumber
            DDG.RGC(ii).Center = DDG.MinCoordinate + (DDG.MaxCoordinate-DDG.MinCoordinate)*rand(DDG.Rng,1,DDG.NumberOfVariables);%Randomly initialize DGCs' center positions inside the boundary
            DDG.RGC(ii).Weight = DDG.MinWeight + (DDG.MaxWeight-DDG.MinWeight)*rand(DDG.Rng);
            switch DDG.Conditioning
                case 0
                    DDG.RGC(ii).Sigma = (DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng))).* ones(1,DDG.NumberOfVariables);
                case 1
                    DDG.RGC(ii).Sigma = DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng,1,DDG.NumberOfVariables));
            end
            switch DDG.Rotation
                case 0
                    DDG.RGC(ii).ThetaMatrix    = zeros(DDG.NumberOfVariables);
                    DDG.RGC(ii).RotationMatrix = eye(DDG.NumberOfVariables);
                case 1
                    DDG.RGC(ii).ThetaMatrix = triu(DDG.MinAngle + (DDG.MaxAngle - DDG.MinAngle) * rand(DDG.Rng,DDG.NumberOfVariables,DDG.NumberOfVariables),1);
                    [DDG.RGC(ii).RotationMatrix] = Rotation(DDG.RGC(ii).ThetaMatrix,DDG.NumberOfVariables);
            end
            DDG.RGC(ii).ShiftSeverity          = DDG.LocalShiftSeverityRange(1)+((DDG.LocalShiftSeverityRange(2)-DDG.LocalShiftSeverityRange(1))* rand(DDG.Rng));%Local shift Severity for relocating the center position of DGC ii
            DDG.RGC(ii).ShiftCorrelationFactor = DDG.RelocationCorrelationRange(1)+((DDG.RelocationCorrelationRange(2)-DDG.RelocationCorrelationRange(1))* rand(DDG.Rng));%Correlation factor for relocating the center position of DGC ii
            tmp                                 = randn(DDG.Rng,1,DDG.NumberOfVariables);
            DDG.RGC(ii).PreviousShiftDirection = tmp/sqrt(sum(tmp.^2)); % Initial shift direction for being used in correlation-based relocation
            DDG.RGC(ii).SigmaSeverity   = DDG.LocalSigmaSeverityRange(1)+((DDG.LocalSigmaSeverityRange(2)-DDG.LocalSigmaSeverityRange(1))* rand(DDG.Rng));%Local sigma Severity for changing the vastness of the basin of attraction of DGC ii
            if DDG.Conditioning==0
                DDG.RGC(ii).SigmaDirection  = ones(1,DDG.NumberOfVariables).*(randi(DDG.Rng,[0, 1], 1, 1) * 2 - 1); % Defines whether ALL sigma values of DGC ii increase or decrease after local changes
            else
                DDG.RGC(ii).SigmaDirection  = randi(DDG.Rng,[0, 1], 1, DDG.NumberOfVariables) * 2 - 1; % Defines whether EACH sigma value of DGC ii increases or decreases after local changes
            end
            DDG.RGC(ii).WeightSeverity  = DDG.LocalWeightSeverityRange(1)+((DDG.LocalWeightSeverityRange(2)-DDG.LocalWeightSeverityRange(1))* rand(DDG.Rng));%Weight Severity for changing the heights of promising regions in the objective space
            DDG.RGC(ii).WeightDirection = randi(DDG.Rng,[0, 1]) * 2 - 1; % Defines whether height increases or decrease after local changes
            DDG.RGC(ii).RotationSeverity  = DDG.LocalRotationSeverityRange(1)+((DDG.LocalRotationSeverityRange(2)-DDG.LocalRotationSeverityRange(1))* rand(DDG.Rng));%Rotation Severity for changing the rotation of the basin of attraction of DGC ii
            DDG.RGC(ii).RotationDirection = triu(randi(DDG.Rng,[0, 1], DDG.NumberOfVariables, DDG.NumberOfVariables) * 2 - 1,1); % Defines whether the rotation for each pair of variables in DGC ii changes clockwise or counter clockwise after local changes
            DDG.RGC(ii).LocalChangeLikelihood      = DDG.LocalTemporalSeverityRange(1)+((DDG.LocalTemporalSeverityRange(2)-DDG.LocalTemporalSeverityRange(1))* rand(DDG.Rng));%Likelihood of change in DGC ii at each function evaluation
            DDG.RGC(ii).DirectionChangeProbabolity = DDG.DirectionChangeProbabolityRange(1)+((DDG.DirectionChangeProbabolityRange(2)-DDG.DirectionChangeProbabolityRange(1))* rand(DDG.Rng));%Likelihood of inverting the direction of changing height, sigma, and angles
        end
    end
    DDG.RGCNumber = UpdatedRGCNumber;
elseif ChangeCode==-2 % Change in the number of variables
    UpdatedVariableNumber = DDG.NumberOfVariables + sign(randn(DDG.Rng))*DDG.VariableNumberChangeSeverity;
    % Boundary check for DGC sigma
    tmp = UpdatedVariableNumber > DDG.MaxNumberOfVariables;
    UpdatedVariableNumber(tmp) = (2*DDG.MaxNumberOfVariables)- UpdatedVariableNumber(tmp);
    tmp = UpdatedVariableNumber < DDG.MinNumberOfVariables;
    UpdatedVariableNumber(tmp) = (2*DDG.MinNumberOfVariables)- UpdatedVariableNumber(tmp);
    if UpdatedVariableNumber<DDG.NumberOfVariables
        VariablesToBeRemoved = randperm(DDG.Rng,DDG.NumberOfVariables,abs(DDG.NumberOfVariables-UpdatedVariableNumber));
        for ii=1:DDG.RGCNumber
            DDG.RGC(ii).Center(VariablesToBeRemoved)=[];
            DDG.RGC(ii).Sigma(VariablesToBeRemoved)=[];
            DDG.RGC(ii).ThetaMatrix(VariablesToBeRemoved, :) = [];
            DDG.RGC(ii).ThetaMatrix(:,VariablesToBeRemoved) = [];
            DDG.RGC(ii).PreviousShiftDirection(VariablesToBeRemoved)=[];
            DDG.RGC(ii).SigmaDirection(VariablesToBeRemoved)=[];
            DDG.RGC(ii).RotationDirection(VariablesToBeRemoved, :) = [];
            DDG.RGC(ii).RotationDirection(:,VariablesToBeRemoved) = [];
            [DDG.RGC(ii).RotationMatrix] = Rotation(DDG.RGC(ii).ThetaMatrix,DDG.NumberOfVariables);
        end
    elseif UpdatedVariableNumber>DDG.NumberOfVariables
        VariablesToBeAdded = sort(randperm(DDG.Rng,UpdatedVariableNumber,abs(DDG.NumberOfVariables-UpdatedVariableNumber)));
        for ii=1:DDG.RGCNumber
            for jj=1:length(VariablesToBeAdded)
                Variable2Add = VariablesToBeAdded(jj);
                % Expand Center and Sigma by shifting elements and inserting new variable
                DDG.RGC(ii).Center = [DDG.RGC(ii).Center(1:Variable2Add-1), (DDG.MinCoordinate + (DDG.MaxCoordinate-DDG.MinCoordinate)*rand(DDG.Rng)), DDG.RGC(ii).Center(Variable2Add:end)];
                tmp = [DDG.RGC(ii).PreviousShiftDirection(1:Variable2Add-1), (randn(DDG.Rng)*mean(DDG.RGC(ii).PreviousShiftDirection)) , DDG.RGC(ii).PreviousShiftDirection(Variable2Add:end)];
                DDG.RGC(ii).PreviousShiftDirection = tmp/sqrt(sum(tmp.^2));
                switch DDG.Conditioning
                    case 0
                        DDG.RGC(ii).Sigma = [DDG.RGC(ii).Sigma, DDG.RGC(ii).Sigma(1)];
                        DDG.RGC(ii).SigmaDirection = [DDG.RGC(ii).SigmaDirection, DDG.RGC(ii).SigmaDirection(1)];
                    case 1
                        DDG.RGC(ii).Sigma = [DDG.RGC(ii).Sigma(1:Variable2Add-1), (DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng))), DDG.RGC(ii).Sigma(Variable2Add:end)];
                        DDG.RGC(ii).SigmaDirection = [DDG.RGC(ii).SigmaDirection(1:Variable2Add-1), randi(DDG.Rng,[0, 1], 1, 1) * 2 - 1, DDG.RGC(ii).SigmaDirection(Variable2Add:end)];
                end
                switch DDG.Rotation
                    case 0
                        DDG.RGC(ii).RotationMatrix = eye(UpdatedVariableNumber);
                        DDG.RGC(ii).ThetaMatrix    = zeros(UpdatedVariableNumber);
                    case 1
                        [row,col] =size(DDG.RGC(ii).ThetaMatrix);
                        DDG.RGC(ii).ThetaMatrix = [DDG.RGC(ii).ThetaMatrix(1:Variable2Add-1, :); (DDG.MinAngle + (DDG.MaxAngle - DDG.MinAngle) * rand(DDG.Rng,1, col)); DDG.RGC(ii).ThetaMatrix(Variable2Add:end, :)];
                        DDG.RGC(ii).ThetaMatrix = [DDG.RGC(ii).ThetaMatrix(:, 1:Variable2Add-1), (DDG.MinAngle + (DDG.MaxAngle - DDG.MinAngle) * rand(DDG.Rng,row+1, 1)), DDG.RGC(ii).ThetaMatrix(:, Variable2Add:end)];
                        DDG.RGC(ii).ThetaMatrix = triu(DDG.RGC(ii).ThetaMatrix,1);
                        DDG.RGC(ii).RotationDirection = [DDG.RGC(ii).RotationDirection(1:Variable2Add-1, :); (randi(DDG.Rng,[0, 1], 1, col) * 2 - 1); DDG.RGC(ii).RotationDirection(Variable2Add:end, :)];
                        DDG.RGC(ii).RotationDirection = [DDG.RGC(ii).RotationDirection(:, 1:Variable2Add-1), (randi(DDG.Rng,[0, 1], row+1,1) * 2 - 1), DDG.RGC(ii).RotationDirection(:, Variable2Add:end)];
                        DDG.RGC(ii).RotationDirection = triu(DDG.RGC(ii).RotationDirection,1);
                end
            end
            [DDG.RGC(ii).RotationMatrix] = Rotation(DDG.RGC(ii).ThetaMatrix,DDG.NumberOfVariables);
        end
    end
    DDG.NumberOfVariables = UpdatedVariableNumber;
elseif ChangeCode==-3 % Change in the number of clusters. Used when the number of clusters is defined by external parts of the system
    UpdatedClusterNumber = DDG.ClusterNumber + sign(randn(DDG.Rng))*DDG.ClusterNumberChangeSeverity;
    % Boundary check for DGC sigma
    tmp = UpdatedClusterNumber > DDG.MaxClusterNumber;
    UpdatedClusterNumber(tmp) = (2*DDG.MaxClusterNumber)- UpdatedClusterNumber(tmp);
    tmp = UpdatedClusterNumber < DDG.MinClusterNumber;
    UpdatedClusterNumber(tmp) = (2*DDG.MinClusterNumber)- UpdatedClusterNumber(tmp);
    DDG.ClusterNumber=UpdatedClusterNumber;
end
%% Updating the MaxWeight value
DDG.MaxWeight         = max(arrayfun(@(x) x.Weight, DDG.RGC));% Value of the largest height
end
%% Rotation matrix generator function
function R = Rotation(teta,Dimension)
R = eye(Dimension);
for p=1 : (Dimension-1)
    for q=(p+1) : (Dimension)
        if teta(p,q)~=0
            G = eye(Dimension);
            G(p,p) = cos(teta(p,q));
            G(q,q) = cos(teta(p,q));
            G(p,q) = -sin(teta(p,q));
            G(q,p) = sin(teta(p,q));
            R = R*G;
        end
    end
end
end