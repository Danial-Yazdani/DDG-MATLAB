%*****Dynamic Dataset Generator (DDG) MATLAB Implementation ver. 1.00******
% Author: X Y
%Last Edited: January 31, 2024
%Title: Main function of DDG
% --------
%Refrence: "Clustering in Dynamic Environments: A Framework for Benchmark
%          Dataset Generation With Heterogeneous Changes"
%
%
% --------
% License:
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: X Y
% e-mail: X DOT Y AT something DOT com
% Copyright notice: (c) 2024 X Y
%**************************************************************************
close all;clear all;clc; %#ok<CLALL>
%% Preparation and initialization
clear DDG;
DDG = DDGinitialization;
rng('shuffle');%Set a random seed for the optimizer
%% Start optimization/clustering method code from here


%Your optimization/clustering method code%

%OfflinePerformance = mean(DDG.CurrentBestSolutionValue);% The performance of an algorithm in a single run.

%% Visualization of data (video) in 2-dimensional space
% v = VideoWriter('DDG_visualization.avi');
% v.FrameRate = 5;
% open(v);
% % The loop
% for ii = 1:25000%Number of function evaluations
%     X = DDG.MinCoordinate + (DDG.MaxCoordinate-DDG.MinCoordinate)*rand(DDG.Rng,1,DDG.ClusterNumber*DDG.NumberOfVariables);
%     [result,DDG] = ClusteringEvaluation(X,DDG);%evaluating a random solution just for activating the environmental changes
%     % Visualize data every N iterations
%     N = 100; % Set N as desired
%     if mod(ii, N) == 0
%         figure('visible', 'off'); % Create an invisible figure
%         scatter(DDG.Data.Dataset(:,2),DDG.Data.Dataset(:,1),10)
%         xlim([-160 160]);
%         ylim([-160 160]);
%         ylabel('x_1');
%         xlabel('x_2');
%         grid on
%         box on
%         set(gcf,'OuterPosition',[150 150 600 550]);
%         legend(['FE = ' num2str(ii) ', DGCs = ' num2str(DDG.DGCNumber)], 'Location', 'northeast');
%         frame = getframe(gcf); % Capture the plot as an image
%         writeVideo(v, frame);
%         close(gcf); % Close the figure
%     end
% end
% close(v);% Close the video writer
%% Visualization of data (image) in 2-dimensional space
% figure;
% scatter(DDG.Data.Dataset(:,2),DDG.Data.Dataset(:,1),10)
% xlim([-160 160]);
% ylim([-160 160]);
% x1lh = ylabel('x_1');
% x2lh = xlabel('x_2');
% grid on
% box on
% set(gcf,'OuterPosition',[150 150 600 550]);
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% x2lh.Position(2) = x2lh.Position(2);
% Figruename = 'Sample1';
% saveas(gcf,Figruename,'epsc')
% saveas(gcf,strcat(Figruename,'.fig'))