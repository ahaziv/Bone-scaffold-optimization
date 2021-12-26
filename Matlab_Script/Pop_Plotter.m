% this script plots the evolving population sampling
clear all; close all; clc;

FilePath = mfilename('fullpath');
FileName = 'Full Data 6-1-2020';
idcs = strfind(FilePath,'\');
ParDir = FilePath(1:idcs(end-1)-1);
addpath(strcat(ParDir,'\Data Files'));
GreedyRawData = xlsread(FileName,strcat('Greedy'));
GeneticRawData = xlsread(FileName,strcat('Genetic'));

%% reading the genetic data
GreedyData = GreedyRawData((GreedyRawData(:,7)==0),1:6);
GreedyData(:,2) = GreedyData(:,2).*GreedyData(:,1);
GreedyData(:,3) = GreedyData(:,3).*GreedyData(:,1);

bounds = GreedyRawData((GreedyRawData(:,7)==2),1:6);
bounds(:,3) = bounds(:,3).* bounds(:,1);
bounds(:,4) = bounds(:,4).* bounds(:,2);
bounds(:,5) = bounds(:,5).* bounds(:,1);
bounds(:,6) = bounds(:,6).* bounds(:,2);

conv_rate = GreedyRawData(1,8);

% cleaning the data
temp_length = 1;
while temp_length ~= length(GreedyData(:,1))
    Avg_Loss = mean(GreedyData(:,6));
    temp_length = length(GreedyData(:,1));
    GreedyData(GreedyData(:,6)>(50*Avg_Loss),:) = [];
end

%% reading the genetic data
GeneticData = GeneticRawData((GeneticRawData(:,7)==0),1:6);
GeneticData = GeneticData((GeneticData(:,5)~=0),1:6);
GeneticData(:,2) = GeneticData(:,2).*GeneticData(:,1);
GeneticData(:,3) = GeneticData(:,3).*GeneticData(:,1);

% cleaning the data
temp_length = 1;
while temp_length ~= length(GeneticData(:,1))
    Avg_Loss = mean(GeneticData(:,6));
    temp_length = length(GeneticData(:,1));
    GeneticData(GeneticData(:,6)>(50*Avg_Loss),:) = [];
end

%% setting the limits for the 3D plots
Limits = zeros(2,3);
Limits(1,:) = 0.9*min([min(GreedyData(:,1:3)); min(GeneticData(:,1:3))]);
Limits(2,:) = 1.1*max([max(GreedyData(:,1:3)); max(GeneticData(:,1:3))]);

%% plotting the 3D results
pos1 = [0.06 0.09 0.43 0.83];
pos2 = [0.55 0.09 0.43 0.83];
pos3 = [0.06 0.09 0.385 0.83];
pos4 = [0.52 0.09 0.48 0.83];
figure
set(gcf, 'Position', [0,0,960,480]);
cmap = parula;

c = -log(GreedyData(:,6)-0.9*min(GreedyData(:,6)));
cn = (c-min(c))/(max(c)-min(c));
cn = ceil(cn*size(cmap,1));
cn = max(cn,1);

subplot ('Position', pos2)
set(gca, 'FontName', 'Times New Roman')
scatter3(GreedyData(:,1),GreedyData(:,2),GreedyData(:,3), 150, cmap(cn(:),:), '.');
title ('Random Search Results', 'fontsize', 16)
caxis([min(GreedyData(:,6)) max(GreedyData(:,6))])
colormap(flipud(cmap))
bar = colorbar;
set(get(bar,'Title'),'String','Loss','fontsize', 12)
grid on

xlabel('FW [\mum]', 'fontsize', 12)
ylabel('FD [\mum]', 'fontsize', 12)
zlabel('ST [\mum]', 'fontsize', 12)
xlim([Limits(1,1) Limits(2,1)])
ylim([Limits(1,2) Limits(2,2)])
zlim([Limits(1,3) Limits(2,3)])
set(gca, 'FontName', 'Times New Roman')

%% a plotting section for genetic data as well
c = -log(GeneticData(:,6)-0.9*min(GeneticData(:,6)));
cn = (c-min(c))/(max(c)-min(c));
cn = ceil(cn*size(cmap,1));
cn = max(cn,1);

subplot ('Position', pos3)
set(gca, 'FontName', 'Times New Roman')
scatter3(GeneticData(:,1),GeneticData(:,2),GeneticData(:,3), 150, cmap(cn(:),:), '.');
title ('Genetic Algorithm Results', 'fontsize', 16)
caxis([min(GeneticData(:,6)) max(GeneticData(:,6))])
colormap(flipud(cmap))
% colorbar
grid on

xlabel('FW [\mum]', 'fontsize', 12)
ylabel('FD [\mum]', 'fontsize', 12)
zlabel('ST [\mum]', 'fontsize', 12)
xlim([Limits(1,1) Limits(2,1)])
ylim([Limits(1,2) Limits(2,2)])
zlim([Limits(1,3) Limits(2,3)])

set(gca, 'FontName', 'Times New Roman')

%% plotting the loss function against r
figure
set(gcf, 'Position', [0,0,960,480]);
subplot ('Position', pos1)
hold on
radius = 0.5*(sqrt((GreedyData(:,2)*10^-6).^2 + (GreedyData(:,3)*10^-6).^2) - (GreedyData(:,1)*10^-6));
scatter(radius, GreedyData(:,6), 40, '.')
radius = 0.5*(sqrt((GeneticData(:,2)*10^-6).^2 + (GeneticData(:,3)*10^-6).^2) - (GeneticData(:,1)*10^-6));
scatter(radius, GeneticData(:,6), 40, '.')

syms r
fplot((-(8.534*10^-8)*log(r-50*10^-6)+(3*r^2)/(2*pi))*10^6)
xlim ([0 700]*10^-6)
ylim ([0.7 1.35])
title ('Loss against Pore Radius', 'fontsize', 16)
xlabel('r [\mum]')
ylabel('Loss rate')
set(gca, 'FontName', 'Times New Roman')
legend('Random Search Loss', 'Genetic Loss', 'Analytic Spherical Loss','fontsize', 12)
legend('Location','northwest')

%% plotting the loss function against r
subplot ('Position', pos2)
hold on
set(gca, 'FontName', 'Times New Roman')
plot(1:length(GreedyData(:,6)), GreedyData(:,6));
plot(1:length(GeneticData(:,6)), GeneticData(:,6));
title ('Algorithm Results', 'fontsize', 16)
grid on

xlim ([0 210])
xlabel('Iteration number', 'fontsize', 12)
ylabel('Loss Rate', 'fontsize', 12)
set(gca, 'FontName', 'Times New Roman')
legend('Random Search Loss', 'Genetic Loss','fontsize', 12)
