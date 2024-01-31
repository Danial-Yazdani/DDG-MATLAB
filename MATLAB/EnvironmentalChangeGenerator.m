%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
% Author: X Y
%Last Edited: January 31, 2024
%Title: Main function of DDG
% --------
%Refrence: "Clustering in Dynamic Environments: A Framework for Benchmark
%          Dataset Generation With Heterogeneous Changes"
%
% Description: This function updates the parameters of the Dynamic Gaussian
% Components (DGCs) to simulate environmental changes within the distributions.
% The parameter 'ChangeCode' specifies the type of change being simulated.
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: X Y
% e-mail: X DOT Y AT something DOT com
% Copyright notice: (c) 2024 X Y
%**************************************************************************
function [DDG] = EnvironmentalChangeGenerator(DDG,ChangeCode)
%The default dynamic used for generating environmental changes is random.
if ChangeCode>0
    %% Gradual local change for the DGC whose ID is equal to ChangeCode
    DGCid = ChangeCode;
    % DGC center update
    tmp = randn(DDG.Rng,1,DDG.NumberOfVariables);
    RandomDirection = tmp/sqrt(sum(tmp.^2));
    SummedVector=(((1-DDG.DGC(DGCid).ShiftCorrelationFactor)*RandomDirection)+ (DDG.DGC(DGCid).ShiftCorrelationFactor*DDG.DGC(DGCid).PreviousShiftDirection));
    RelocationDirection = SummedVector/ sqrt(sum(SummedVector.^2));% Normalized to create a unit vector
    UpdatedDGCsPosition  = DDG.DGC(DGCid).Center + (RelocationDirection*(abs(randn(DDG.Rng))*DDG.DGC(DGCid).ShiftSeverity));
    % Boundary check for DGC center
    tmp = UpdatedDGCsPosition > DDG.MaxCoordinate;
    UpdatedDGCsPosition(tmp) = (2*DDG.MaxCoordinate)- UpdatedDGCsPosition(tmp);
    tmp = UpdatedDGCsPosition < DDG.MinCoordinate;
    UpdatedDGCsPosition(tmp) = (2*DDG.MinCoordinate)- UpdatedDGCsPosition(tmp);
    RelocationVectore = UpdatedDGCsPosition-DDG.DGC(DGCid).Center;
    DDG.DGC(DGCid).PreviousShiftDirection = (RelocationVectore)./ sqrt(sum(RelocationVectore.^2));%Storing the last relocation direction
    DDG.DGC(DGCid).Center = UpdatedDGCsPosition;
    % Weight update
    if rand(DDG.Rng)<DDG.DGC(DGCid).DirectionChangeProbabolity
        DDG.DGC(DGCid).WeightDirection = -DDG.DGC(DGCid).WeightDirection;%Update heigh change direction based on the probability
    end
    UpdatedWeight = DDG.DGC(DGCid).Weight + (abs(randn(DDG.Rng))*DDG.DGC(DGCid).WeightSeverity*DDG.DGC(DGCid).WeightDirection);
    % Boundary check for DGC height
    tmp = UpdatedWeight > DDG.MaxWeight;
    DDG.DGC(DGCid).WeightDirection(tmp) = -DDG.DGC(DGCid).WeightDirection(tmp); % Changing direction if reached the boundary
    UpdatedWeight(tmp) = (2*DDG.MaxWeight)- UpdatedWeight(tmp);
    tmp = UpdatedWeight < DDG.MinWeight;
    DDG.DGC(DGCid).WeightDirection(tmp) = -DDG.DGC(DGCid).WeightDirection(tmp); % Changing direction if reached the boundary
    UpdatedWeight(tmp) = (2*DDG.MinWeight)- UpdatedWeight(tmp);
    DDG.DGC(DGCid).Weight = UpdatedWeight;
    % Sigma update
    if DDG.Conditioning == 0%Keeps condition number at 1
        if rand(DDG.Rng)<DDG.DGC(DGCid).DirectionChangeProbabolity
            DDG.DGC(DGCid).SigmaDirection = -DDG.DGC(DGCid).SigmaDirection;
        end
        UpdatedSigma = DDG.DGC(DGCid).Sigma + ( (ones(1,DDG.NumberOfVariables).*abs(randn(DDG.Rng))) .* DDG.DGC(DGCid).SigmaSeverity .* DDG.DGC(DGCid).SigmaDirection);
    else
        invertFlags = rand(DDG.Rng,1,DDG.NumberOfVariables) < DDG.DGC(DGCid).DirectionChangeProbabolity;%Update the elements of sigma change direction based on the probability
        DDG.DGC(DGCid).SigmaDirection(invertFlags) = -DDG.DGC(DGCid).SigmaDirection(invertFlags);
        UpdatedSigma = DDG.DGC(DGCid).Sigma + (abs(randn(DDG.Rng,1,DDG.NumberOfVariables))*DDG.DGC(DGCid).SigmaSeverity.*DDG.DGC(DGCid).SigmaDirection);
    end
    % Boundary check for DGC sigma
    tmp = UpdatedSigma > DDG.MaxSigma;
    DDG.DGC(DGCid).SigmaDirection(tmp) = -DDG.DGC(DGCid).SigmaDirection(tmp); % Changing direction if reached the boundary
    UpdatedSigma(tmp) = (2*DDG.MaxSigma)- UpdatedSigma(tmp);
    tmp = UpdatedSigma < DDG.MinSigma;
    DDG.DGC(DGCid).SigmaDirection(tmp) = -DDG.DGC(DGCid).SigmaDirection(tmp); % Changing direction if reached the boundary
    UpdatedSigma(tmp) = (2*DDG.MinSigma)- UpdatedSigma(tmp);
    DDG.DGC(DGCid).Sigma = UpdatedSigma;
    % Angle update
    if DDG.Rotation==1
        invertFlags = triu((rand(DDG.Rng,DDG.NumberOfVariables,DDG.NumberOfVariables) < DDG.DGC(DGCid).DirectionChangeProbabolity),1);
        DDG.DGC(DGCid).RotationDirection(invertFlags) = -DDG.DGC(DGCid).RotationDirection(invertFlags);
        UpdatedAngles = DDG.DGC(DGCid).ThetaMatrix + (triu(abs(randn(DDG.Rng,DDG.NumberOfVariables,DDG.NumberOfVariables)),1).*DDG.DGC(DGCid).RotationSeverity.*DDG.DGC(DGCid).RotationDirection);
        % Boundary check for angles
        tmp = UpdatedAngles > DDG.MaxAngle;
        DDG.DGC(DGCid).RotationDirection(tmp) = -DDG.DGC(DGCid).RotationDirection(tmp); % Changing direction if reached the boundary
        UpdatedAngles(tmp) = (2*DDG.MaxAngle)- UpdatedAngles(tmp);
        tmp = UpdatedAngles < DDG.MinAngle;
        DDG.DGC(DGCid).RotationDirection(tmp) = -DDG.DGC(DGCid).RotationDirection(tmp); % Changing direction if reached the boundary
        UpdatedAngles(tmp) = (2*DDG.MinAngle)- UpdatedAngles(tmp);
        DDG.DGC(DGCid).ThetaMatrix = UpdatedAngles;
        [DDG.DGC(DGCid).RotationMatrix] = Rotation(DDG.DGC(DGCid).ThetaMatrix,DDG.NumberOfVariables);
    end
elseif ChangeCode==0
    %% Severe global change for all DGCs
    RandStream.setGlobalStream(DDG.Rng); % Set the global RandStream to DDG.Rng since betarand does not accpet DDG.Rng as RandStream.
    for DGCid=1 : DDG.DGCNumber
        % DGC center update
        tmp = randn(DDG.Rng,1,DDG.NumberOfVariables);
        RandomDirection = tmp/sqrt(sum(tmp.^2));
        UpdatedDGCsPosition  = DDG.DGC(DGCid).Center + (RandomDirection*DDG.GlobalShiftSeverityValue*(2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl)-1));
        % Boundary check for DGC center
        tmp = UpdatedDGCsPosition > DDG.MaxCoordinate;
        UpdatedDGCsPosition(tmp) = (2*DDG.MaxCoordinate)- UpdatedDGCsPosition(tmp);
        tmp = UpdatedDGCsPosition < DDG.MinCoordinate;
        UpdatedDGCsPosition(tmp) = (2*DDG.MinCoordinate)- UpdatedDGCsPosition(tmp);
        DDG.DGC(DGCid).Center = UpdatedDGCsPosition;
        % Weight update
        UpdatedWeight = DDG.DGC(DGCid).Weight + (DDG.GlobalWeightSeverityValue*(2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl)-1));
        % Boundary check for DGC height
        tmp = UpdatedWeight > DDG.MaxWeight;
        UpdatedWeight(tmp) = (2*DDG.MaxWeight)- UpdatedWeight(tmp);
        tmp = UpdatedWeight < DDG.MinWeight;
        UpdatedWeight(tmp) = (2*DDG.MinWeight)- UpdatedWeight(tmp);
        DDG.DGC(DGCid).Weight = UpdatedWeight;
        % Sigma update
        if DDG.Conditioning == 0%Keeps condition number at 1
            UpdatedSigma = DDG.DGC(DGCid).Sigma + ( (ones(1,DDG.NumberOfVariables).* (2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl)-1)) .* DDG.GlobalSigmaSeverityValue);
        else
            UpdatedSigma = DDG.DGC(DGCid).Sigma + ( (2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl,1,DDG.NumberOfVariables)-1) .* DDG.GlobalSigmaSeverityValue);
        end
        % Boundary check for DGC sigma
        tmp = UpdatedSigma > DDG.MaxSigma;
        UpdatedSigma(tmp) = (2*DDG.MaxSigma)- UpdatedSigma(tmp);
        tmp = UpdatedSigma < DDG.MinSigma;
        UpdatedSigma(tmp) = (2*DDG.MinSigma)- UpdatedSigma(tmp);
        DDG.DGC(DGCid).Sigma = UpdatedSigma;
        % Angle update
        if DDG.Rotation==1
            UpdatedAngles = DDG.DGC(DGCid).ThetaMatrix + (DDG.GlobalAngleSeverityValue*triu((2*betarnd(DDG.GlobalSeverityControl,DDG.GlobalSeverityControl,DDG.NumberOfVariables,DDG.NumberOfVariables)-1),1));
            % Boundary check for angles
            tmp = UpdatedAngles > DDG.MaxAngle;
            UpdatedAngles(tmp) = (2*DDG.MaxAngle)- UpdatedAngles(tmp);
            tmp = UpdatedAngles < DDG.MinAngle;
            UpdatedAngles(tmp) = (2*DDG.MinAngle)- UpdatedAngles(tmp);
            DDG.DGC(DGCid).ThetaMatrix = UpdatedAngles;
            [DDG.DGC(DGCid).RotationMatrix] = Rotation(DDG.DGC(DGCid).ThetaMatrix,DDG.NumberOfVariables);
        end
    end
    rng('shuffle');%Shuffle the global randstream
elseif ChangeCode==-1
    %% Change in the number of DGCs
    UpdatedDGCNumber = DDG.DGCNumber + (randi(DDG.Rng,[0, 1]) * 2 - 1)*DDG.DGCNumberChangeSeverity;
    % Boundary check for DGC sigma
    tmp = UpdatedDGCNumber > DDG.MaxDGCNumber;
    UpdatedDGCNumber(tmp) = (2*DDG.MaxDGCNumber)- UpdatedDGCNumber(tmp);
    tmp = UpdatedDGCNumber < DDG.MinDGCNumber;
    UpdatedDGCNumber(tmp) = (2*DDG.MinDGCNumber)- UpdatedDGCNumber(tmp);
    if UpdatedDGCNumber<DDG.DGCNumber
        DGCsToBeRemoved = randperm(DDG.Rng,DDG.DGCNumber,abs(DDG.DGCNumber-UpdatedDGCNumber));
        DDG.DGC(DGCsToBeRemoved)=[];% Removing randomly chosen DGCs
    elseif UpdatedDGCNumber>DDG.DGCNumber
        for ii=DDG.DGCNumber+1 : UpdatedDGCNumber
            DDG.DGC(ii).Center = DDG.MinCoordinate + (DDG.MaxCoordinate-DDG.MinCoordinate)*rand(DDG.Rng,1,DDG.NumberOfVariables);%Randomly initialize DGCs' center positions inside the boundary
            DDG.DGC(ii).Weight = DDG.MinWeight + (DDG.MaxWeight-DDG.MinWeight)*rand(DDG.Rng);
            switch DDG.Conditioning
                case 0
                    DDG.DGC(ii).Sigma = (DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng))).* ones(1,DDG.NumberOfVariables);
                case 1
                    DDG.DGC(ii).Sigma = DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng,1,DDG.NumberOfVariables));
            end
            switch DDG.Rotation
                case 0
                    DDG.DGC(ii).ThetaMatrix    = zeros(DDG.NumberOfVariables);
                    DDG.DGC(ii).RotationMatrix = eye(DDG.NumberOfVariables);
                case 1
                    DDG.DGC(ii).ThetaMatrix = triu(DDG.MinAngle + (DDG.MaxAngle - DDG.MinAngle) * rand(DDG.Rng,DDG.NumberOfVariables,DDG.NumberOfVariables),1);
                    [DDG.DGC(ii).RotationMatrix] = Rotation(DDG.DGC(ii).ThetaMatrix,DDG.NumberOfVariables);
            end
            DDG.DGC(ii).ShiftSeverity          = DDG.LocalShiftSeverityRange(1)+((DDG.LocalShiftSeverityRange(2)-DDG.LocalShiftSeverityRange(1))* rand(DDG.Rng));%Local shift Severity for relocating the center position of DGC ii
            DDG.DGC(ii).ShiftCorrelationFactor = DDG.RelocationCorrelationRange(1)+((DDG.RelocationCorrelationRange(2)-DDG.RelocationCorrelationRange(1))* rand(DDG.Rng));%Correlation factor for relocating the center position of DGC ii
            tmp                                 = randn(DDG.Rng,1,DDG.NumberOfVariables);
            DDG.DGC(ii).PreviousShiftDirection = tmp/sqrt(sum(tmp.^2)); % Initial shift direction for being used in correlation-based relocation
            DDG.DGC(ii).SigmaSeverity   = DDG.LocalSigmaSeverityRange(1)+((DDG.LocalSigmaSeverityRange(2)-DDG.LocalSigmaSeverityRange(1))* rand(DDG.Rng));%Local sigma Severity for changing the vastness of the basin of attraction of DGC ii
            if DDG.Conditioning==0
                DDG.DGC(ii).SigmaDirection  = ones(1,DDG.NumberOfVariables).*(randi(DDG.Rng,[0, 1], 1, 1) * 2 - 1); % Defines whether ALL sigma values of DGC ii increase or decrease after local changes
            else
                DDG.DGC(ii).SigmaDirection  = randi(DDG.Rng,[0, 1], 1, DDG.NumberOfVariables) * 2 - 1; % Defines whether EACH sigma value of DGC ii increases or decreases after local changes
            end
            DDG.DGC(ii).WeightSeverity  = DDG.LocalWeightSeverityRange(1)+((DDG.LocalWeightSeverityRange(2)-DDG.LocalWeightSeverityRange(1))* rand(DDG.Rng));%Weight Severity for changing the heights of promising regions in the objective space
            DDG.DGC(ii).WeightDirection = randi(DDG.Rng,[0, 1]) * 2 - 1; % Defines whether height increases or decrease after local changes
            DDG.DGC(ii).RotationSeverity  = DDG.LocalRotationSeverityRange(1)+((DDG.LocalRotationSeverityRange(2)-DDG.LocalRotationSeverityRange(1))* rand(DDG.Rng));%Rotation Severity for changing the rotation of the basin of attraction of DGC ii
            DDG.DGC(ii).RotationDirection = triu(randi(DDG.Rng,[0, 1], DDG.NumberOfVariables, DDG.NumberOfVariables) * 2 - 1,1); % Defines whether the rotation for each pair of variables in DGC ii changes clockwise or counter clockwise after local changes
            DDG.DGC(ii).LocalChangeLikelihood      = DDG.LocalTemporalSeverityRange(1)+((DDG.LocalTemporalSeverityRange(2)-DDG.LocalTemporalSeverityRange(1))* rand(DDG.Rng));%Likelihood of change in DGC ii at each function evaluation
            DDG.DGC(ii).DirectionChangeProbabolity = DDG.DirectionChangeProbabolityRange(1)+((DDG.DirectionChangeProbabolityRange(2)-DDG.DirectionChangeProbabolityRange(1))* rand(DDG.Rng));%Likelihood of inverting the direction of changing height, sigma, and angles
        end
    end
    DDG.DGCNumber = UpdatedDGCNumber;
elseif ChangeCode==-2
    %% Change in the number of variables
    UpdatedVariableNumber = DDG.NumberOfVariables + (randi(DDG.Rng,[0, 1]) * 2 - 1)*DDG.VariableNumberChangeSeverity;
    % Boundary check for DGC sigma
    tmp = UpdatedVariableNumber > DDG.MaxNumberOfVariables;
    UpdatedVariableNumber(tmp) = (2*DDG.MaxNumberOfVariables)- UpdatedVariableNumber(tmp);
    tmp = UpdatedVariableNumber < DDG.MinNumberOfVariables;
    UpdatedVariableNumber(tmp) = (2*DDG.MinNumberOfVariables)- UpdatedVariableNumber(tmp);
    if UpdatedVariableNumber<DDG.NumberOfVariables
        VariablesToBeRemoved = randperm(DDG.Rng,DDG.NumberOfVariables,abs(DDG.NumberOfVariables-UpdatedVariableNumber));
        for ii=1:DDG.DGCNumber
            DDG.DGC(ii).Center(VariablesToBeRemoved)=[];
            DDG.DGC(ii).Sigma(VariablesToBeRemoved)=[];
            DDG.DGC(ii).ThetaMatrix(VariablesToBeRemoved, :) = [];
            DDG.DGC(ii).ThetaMatrix(:,VariablesToBeRemoved) = [];
            DDG.DGC(ii).PreviousShiftDirection(VariablesToBeRemoved)=[];
            DDG.DGC(ii).SigmaDirection(VariablesToBeRemoved)=[];
            DDG.DGC(ii).RotationDirection(VariablesToBeRemoved, :) = [];
            DDG.DGC(ii).RotationDirection(:,VariablesToBeRemoved) = [];
            [DDG.DGC(ii).RotationMatrix] = Rotation(DDG.DGC(ii).ThetaMatrix,DDG.NumberOfVariables);
        end
    elseif UpdatedVariableNumber>DDG.NumberOfVariables
        VariablesToBeAdded = sort(randperm(DDG.Rng,UpdatedVariableNumber,abs(DDG.NumberOfVariables-UpdatedVariableNumber)));
        for ii=1:DDG.DGCNumber
            for jj=1:length(VariablesToBeAdded)
                Variable2Add = VariablesToBeAdded(jj);
                % Expand Center and Sigma by shifting elements and inserting new variable
                DDG.DGC(ii).Center = [DDG.DGC(ii).Center(1:Variable2Add-1), (DDG.MinCoordinate + (DDG.MaxCoordinate-DDG.MinCoordinate)*rand(DDG.Rng)), DDG.DGC(ii).Center(Variable2Add:end)];
                tmp = [DDG.DGC(ii).PreviousShiftDirection(1:Variable2Add-1), (randn(DDG.Rng)*mean(DDG.DGC(ii).PreviousShiftDirection)) , DDG.DGC(ii).PreviousShiftDirection(Variable2Add:end)];
                DDG.DGC(ii).PreviousShiftDirection = tmp/sqrt(sum(tmp.^2));
                switch DDG.Conditioning
                    case 0
                        DDG.DGC(ii).Sigma = [DDG.DGC(ii).Sigma, DDG.DGC(ii).Sigma(1)];
                        DDG.DGC(ii).SigmaDirection = [DDG.DGC(ii).SigmaDirection, DDG.DGC(ii).SigmaDirection(1)];
                    case 1
                        DDG.DGC(ii).Sigma = [DDG.DGC(ii).Sigma(1:Variable2Add-1), (DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng))), DDG.DGC(ii).Sigma(Variable2Add:end)];
                        DDG.DGC(ii).SigmaDirection = [DDG.DGC(ii).SigmaDirection(1:Variable2Add-1), randi(DDG.Rng,[0, 1], 1, 1) * 2 - 1, DDG.DGC(ii).SigmaDirection(Variable2Add:end)];
                end
                switch DDG.Rotation
                    case 0
                        DDG.DGC(ii).RotationMatrix = eye(UpdatedVariableNumber);
                        DDG.DGC(ii).ThetaMatrix    = zeros(UpdatedVariableNumber);
                    case 1
                        [row,col] =size(DDG.DGC(ii).ThetaMatrix);
                        DDG.DGC(ii).ThetaMatrix = [DDG.DGC(ii).ThetaMatrix(1:Variable2Add-1, :); (DDG.MinAngle + (DDG.MaxAngle - DDG.MinAngle) * rand(DDG.Rng,1, col)); DDG.DGC(ii).ThetaMatrix(Variable2Add:end, :)];
                        DDG.DGC(ii).ThetaMatrix = [DDG.DGC(ii).ThetaMatrix(:, 1:Variable2Add-1), (DDG.MinAngle + (DDG.MaxAngle - DDG.MinAngle) * rand(DDG.Rng,row+1, 1)), DDG.DGC(ii).ThetaMatrix(:, Variable2Add:end)];
                        DDG.DGC(ii).ThetaMatrix = triu(DDG.DGC(ii).ThetaMatrix,1);
                        DDG.DGC(ii).RotationDirection = [DDG.DGC(ii).RotationDirection(1:Variable2Add-1, :); (randi(DDG.Rng,[0, 1], 1, col) * 2 - 1); DDG.DGC(ii).RotationDirection(Variable2Add:end, :)];
                        DDG.DGC(ii).RotationDirection = [DDG.DGC(ii).RotationDirection(:, 1:Variable2Add-1), (randi(DDG.Rng,[0, 1], row+1,1) * 2 - 1), DDG.DGC(ii).RotationDirection(:, Variable2Add:end)];
                        DDG.DGC(ii).RotationDirection = triu(DDG.DGC(ii).RotationDirection,1);
                end
            end
            [DDG.DGC(ii).RotationMatrix] = Rotation(DDG.DGC(ii).ThetaMatrix,DDG.NumberOfVariables);
        end
    end
    DDG.NumberOfVariables = UpdatedVariableNumber;
elseif ChangeCode==-3
    %% Change in the number of clusters. Used when the number of clusters is defined by external parts of the system
    UpdatedClusterNumber = DDG.ClusterNumber + (randi(DDG.Rng,[0, 1]) * 2 - 1)*DDG.ClusterNumberChangeSeverity;
    % Boundary check for DGC sigma
    tmp = UpdatedClusterNumber > DDG.MaxClusterNumber;
    UpdatedClusterNumber(tmp) = (2*DDG.MaxClusterNumber)- UpdatedClusterNumber(tmp);
    tmp = UpdatedClusterNumber < DDG.MinClusterNumber;
    UpdatedClusterNumber(tmp) = (2*DDG.MinClusterNumber)- UpdatedClusterNumber(tmp);
    DDG.ClusterNumber=UpdatedClusterNumber;
end
end
%% Generating rotation matrix for a DGC based on matrix Theta
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