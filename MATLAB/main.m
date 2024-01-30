%**************Dynamic Dataset Generator (DDG) ver. 1.00*******************
%Author: Danial Yazdani
%Last Edited: January 17, 2024
%Title: Main function of DDG
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
% Copyright notice: (c) 2024 Danial Yazdani
%**************************************************************************
% close all;clear all;clc; %#ok<CLALL>
%% Preparation and initialization
clear DDG;
DDG = DDGinitialization;
rng('shuffle');%Set a random seed for the optimizer
% for ii=1:1000
%     X = DDG.MinCoordinate + (DDG.MaxCoordinate-DDG.MinCoordinate)*rand(DDG.Rng,1,DDG.ClusterNumber*DDG.NumberOfVariables);
%     [result,DDG] = ClusteringEvaluation(X,DDG);
% end
%% Visualization
T = DDG.MinCoordinate : ( DDG.MaxCoordinate-DDG.MinCoordinate)/200 :  DDG.MaxCoordinate;
L=length(T);
F=zeros(L);
for i=1:L
    for j=1:L
        [F(i,j),DDG] = DynamicLandscapeFunction([T(i), T(j)],DDG);
    end
end
% surf plot
figure;
f = surf(T,T,F,'FaceAlpha',0.8, 'EdgeColor', 'none');
x1lh = ylabel('x_1');
x2lh = xlabel('x_2');
zlabel('f(x_1,x_2)')
colormap jet
grid on
set(gcf,'OuterPosition',[150 150 600 550]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
x2lh.Position(2) = x2lh.Position(2);
Figruename = 'GMPBsurf1';
saveas(gcf,Figruename,'epsc')
saveas(gcf,strcat(Figruename,'.fig'))
% contour plot
figure;
contour(T,T,F,100,'ZLocation','zmin');
x1lh = ylabel('x_1');
x2lh = xlabel('x_2');
zlabel('f(x_1,x_2)')
colormap jet
grid on
set(gcf,'OuterPosition',[150 150 600 550]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
x2lh.Position(2) = x2lh.Position(2);
Figruename = 'GMPBcontour1';
saveas(gcf,Figruename,'epsc')
saveas(gcf,strcat(Figruename,'.fig'))
% Sampling process/Data generation
% switch DataGenerationApproach
%     case 1
%         Data.list = PDFsampling(Data.SampleSize,DDG,F);
%     case 2
DDG = DFRSsampling(DDG.Data.SampleSize,DDG);
% end
figure;
scatter(DDG.Data.Dataset(:,2),DDG.Data.Dataset(:,1),10)
xlim([-100 100]);
ylim([-100 100]);
x1lh = ylabel('x_1');
x2lh = xlabel('x_2');
grid on
box on
set(gcf,'OuterPosition',[150 150 600 550]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
x2lh.Position(2) = x2lh.Position(2);
Figruename = 'Sample1';
saveas(gcf,Figruename,'epsc')
saveas(gcf,strcat(Figruename,'.fig'))