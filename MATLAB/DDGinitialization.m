%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
% Author: X Y
%Last Edited: January 31, 2024
%Title: Main function of DDG
% --------
%Refrence: "Clustering in Dynamic Environments: A Framework for Benchmark 
%          Dataset Generation With Heterogeneous Changes"
%
% --------
% Description: This function initializes the Dynamic Dataset Generator (DDG), 
% setting the stage for generating diverse dynamic scenarios. Users can 
% customize the DDG's behavior by adjusting its parameters directly within 
% this function. This flexibility allows for the tailored creation of 
% scenarios to meet specific research needs or objectives.
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: X Y
% e-mail: X DOT Y AT something DOT com
% Copyright notice: (c) 2024 X Y
%************************************************************************** 
function [DDG] = DDGinitialization
DDG          = [];
DDG.Seed     = 321; 
DDG.Rng      = RandStream('mt19937ar', 'Seed', DDG.Seed); % DDG uses random numbers to generate problem instances which affect some problem characteristics. An identical random seed should be used for all runs to generate identical problem instances for each experiment.
DDG.MaxEvals = 500000; % Maximum function evaluation number
%% Number of DGCs, variables, and clusters
DDG.MinNumberOfVariables = 2;
DDG.MaxNumberOfVariables = 5;
DDG.NumberOfVariables    = 2;% Activate the next line if you want to randomly initialize the number of variables
%DDG.NumberOfVariables    = randi(DDG.Rng,[DDG.MinNumberOfVariables, DDG.MaxNumberOfVariables]); % Initial number of variables in the dataset
DDG.MinDGCNumber         = 2; 
DDG.MaxDGCNumber         = 10; 
DDG.DGCNumber            = 5;% Activate the next line if you want to randomly initialize the number of DGCs
% DDG.DGCNumber            = randi(DDG.Rng,[DDG.MinDGCNumber, DDG.MaxDGCNumber]);% Initial number of DGCs in the landscape
DDG.MinClusterNumber     = 2; 
DDG.MaxClusterNumber     = 10; 
DDG.ClusterNumber        = 5;% Activate the next line if you want to randomly initialize the number of clusters
% DDG.ClusterNumber        = randi(DDG.Rng,[DDG.MinClusterNumber, DDG.MaxClusterNumber]);% Initial number of clusters.
%% Initializing the center positions of DGCs. 
DDG.MinCoordinate     = -100; % Used for bounding the center (mean) positions of DGCs
DDG.MaxCoordinate     = 100; % Used for bounding the center (mean) positions of DGCs
for ii=1:DDG.DGCNumber
    DDG.DGC(ii).Center = DDG.MinCoordinate + (DDG.MaxCoordinate-DDG.MinCoordinate)*rand(DDG.Rng,1,DDG.NumberOfVariables);%Randomly initialize DGCs' center positions inside the boundary
end
%% Defining the weight values of the DGCs
DDG.MinWeight = 1;
DDG.MaxWeight = 3;
for ii=1:DDG.DGCNumber
    DDG.DGC(ii).Weight = DDG.MinWeight + (DDG.MaxWeight-DDG.MinWeight)*rand(DDG.Rng);
end
%% Defining sigma (standard deviation) of DGCs (Configuring condition number)
DDG.MinSigma = 7;
DDG.MaxSigma = 20;
DDG.Conditioning  = 1;% (0) Condition number is 1 for all DGCs but the sigma values are different from a DGC to another.
                      % (1) Condition number is random for all DGCs.
switch DDG.Conditioning
    case 0
        for ii=1:DDG.DGCNumber
            DDG.DGC(ii).Sigma = (DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng))).* ones(1,DDG.NumberOfVariables);
        end
    case 1
        for ii=1:DDG.DGCNumber
            DDG.DGC(ii).Sigma = DDG.MinSigma + ((DDG.MaxSigma-DDG.MinSigma)*rand(DDG.Rng,1,DDG.NumberOfVariables));
        end
    otherwise
        warning('Wrong number is chosen for DDG.Conditioning.')
end
%% Angles and rotatios
DDG.MinAngle = -pi;
DDG.MaxAngle = pi;
DDG.Rotation = 1;%(0) Without rotation
                 %(1) Random Rotation for all DGCs==> Rotation with random angles for each plane for each DGC
switch DDG.Rotation
    case 0
        DDG.RotationMatrix = NaN(DDG.NumberOfVariables,DDG.NumberOfVariables,DDG.DGCNumber);
        DDG.ThetaMatrix = NaN(DDG.NumberOfVariables,DDG.NumberOfVariables,DDG.DGCNumber);
        for ii=1 : DDG.DGCNumber
            DDG.DGC(ii).RotationMatrix = eye(DDG.NumberOfVariables);
            DDG.DGC(ii).ThetaMatrix    = zeros(DDG.NumberOfVariables);
        end
    case 1
        for ii=1 : DDG.DGCNumber
            DDG.ThetaMatrix = zeros(DDG.NumberOfVariables);
            ThetaMatrix = zeros(DDG.NumberOfVariables);
            upperTriangle = triu(true(DDG.NumberOfVariables), 1);
            ThetaMatrix(upperTriangle) = DDG.MinAngle + (DDG.MaxAngle - DDG.MinAngle) * rand(DDG.Rng,sum(upperTriangle(:)), 1);
            DDG.DGC(ii).ThetaMatrix = ThetaMatrix;
            [DDG.DGC(ii).RotationMatrix] = Rotation(DDG.DGC(ii).ThetaMatrix,DDG.NumberOfVariables);
        end
    otherwise
        warning('Wrong number is chosen for DDG.Rotation.')
end
%% Change severity values for Gradual local changes for each DGC
%For parameters that are not going to be impacted in environmental changes (i.e., remain fixed over time), set the severity range to [0,0].
DDG.LocalShiftSeverityRange        = [0.1,0.2];
DDG.RelocationCorrelationRange     = [0.99,0.995];
DDG.LocalSigmaSeverityRange        = [0.05,0.1];
DDG.LocalWeightSeverityRange       = [0.02,0.05];
DDG.LocalRotationSeverityRange     = [pi/360,pi/180];
DDG.DirectionChangeProbabolityRange= [0.02,0.05];
DDG.LocalTemporalSeverityRange     = [0.05,0.1];
for ii=1 : DDG.DGCNumber
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
%% Change severity values for severe changes in the parameters of all DGCs
DDG.GlobalShiftSeverityValue  = 5;
DDG.GlobalSigmaSeverityValue  = 5;
DDG.GlobalWeightSeverityValue = 0.5;
DDG.GlobalAngleSeverityValue  = pi/4;
DDG.GlobalSeverityControl     = 0.1;% The values of alpha and beta in Beta-distribution used in global changes. The values must be 0<alpha=beta<=1, and smaller values result in more heavy-tail distributions.
DDG.GlobalChangeLikelihood    = 0.0008;
%% Parameters for changing the numbers of variables, DGCs, and cluster centers
DDG.DGCNumberChangeSeverity        = 1;
DDG.VariableNumberChangeSeverity   = 1;
DDG.ClusterNumberChangeSeverity    = 1;
DDG.DGCNumberChangeLikelihood      = 0.0004;
DDG.VariableNumberChangeLikelihood = 0;
DDG.ClusterNumberChangeLikelihood  = 0;
%% Parameters used for storing the results
DDG.BestValueAtEachFE        = inf(1,DDG.MaxEvals);
DDG.FE                       = 0;
DDG.CurrentBestSolution      = [];
DDG.CurrentBestSolutionValue = inf; % Set it to inf for minimization problems and -inf for maximization
%% Defining dataset and sampling parameters
DDG.Data.Size      = 1000;
DDG.Data.FrequentSamplingLikelihood = 0.1;% The likelihood of Incremental Sampling
DDG.Data.IncrementalSamplingSize = ceil(DDG.Data.Size*0.05);% Define the percentage of dataset to be replaced by new samples
DDG.Data.Dataset = NaN(DDG.Data.Size,DDG.NumberOfVariables);
DDG = DataGeneration(DDG.Data.Size,DDG);
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