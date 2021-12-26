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
winner = GreedyRawData((GreedyRawData(:,7)==1),1:3);
winner(1,:) = []; 
winner(:,2) = winner(:,2).*winner(:,1);
winner(:,3) = winner(:,3).*winner(:,1);

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
GeneticData = GeneticRawData((GeneticRawData(:,7)==1),1:6);
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

%% plotting FW-FD
figure
set(gcf, 'Position', [0,0,960,480]);
pos1 = [0.06 0.09 0.43 0.83];
subplot ('Position', pos1)
hold on

scatter(GreedyData(:,1), GreedyData(:,2),150,'.')
p1 = fplot(@(FW) 1.1*FW, '--','Color',[0.8350, 0.0780, 0.0840],'Linewidth',1.5);
fplot(@(FW) 3*FW, '--','Color',[0.8350, 0.0780, 0.0840],'Linewidth',1.5)
for ii = 1:length(bounds(:,1))
    four_pts = [bounds(ii,1) bounds(ii,3)
                bounds(ii,1) bounds(ii,4)
                bounds(ii,2) bounds(ii,4)
                bounds(ii,2) bounds(ii,3)
                bounds(ii,1) bounds(ii,3)];    
    
    p2 = plot(four_pts(:,1), four_pts(:,2), 'Color', [0, 0.6, 0],'Linewidth',0.75);
    temp = sprintf("%d", ii);
%     text(four_pts(3,1), four_pts(2,2),temp, 'FontName', 'Times New Roman', 'fontsize', 12)
    text(mean(four_pts(2:3,1)), four_pts(2,2),temp, 'FontName', 'Times New Roman', 'fontsize', 12)
    
end
scatter(winner(:,1), winner(:,2), 200, '.', 'MarkerEdgeColor',[0.9290, 0.6940, 0.1250])
title ('Search Boundaries FW-FD axis', 'fontsize', 16)
legend([p1,p2],'Feasibility Boundary','Search Areas', 'fontsize', 12)
legend boxoff
xlim([0.9*bounds(1,1) 1.05*bounds(1,2)])
ylim([0.9*bounds(1,3) 1.05*bounds(1,4)])
xlabel('FW [\mum]', 'fontsize', 12)
ylabel('FD [\mum]', 'fontsize', 12)
set(gca, 'FontName', 'Times New Roman')

%% plotting FW-ST
pos2 = [0.55 0.09 0.43 0.83];
subplot ('Position', pos2)
hold on

scatter(GreedyData(:,1), GreedyData(:,3),150,'.')
p1 = fplot(@(FW) 0.9*FW, '--','Color',[0.8350, 0.0780, 0.0840],'Linewidth',1.5);
fplot(@(FW) 0.6*FW, '--','Color',[0.8350, 0.0780, 0.0840],'Linewidth',1.5)
for ii = 1:length(bounds(:,1))
    four_pts = [bounds(ii,1) bounds(ii,5)
                bounds(ii,1) bounds(ii,6)
                bounds(ii,2) bounds(ii,6)
                bounds(ii,2) bounds(ii,5)
                bounds(ii,1) bounds(ii,5)];
    
    p2 = plot(four_pts(:,1), four_pts(:,2),'Color', [0, 0.6, 0],'Linewidth',0.75);
    temp = sprintf("%d", ii);        
%     if mod(ii,2)
%         text(1*four_pts(3,1), four_pts(2,2),temp, 'FontName', 'Times New Roman', 'fontsize', 10)
%     else
%         text(1*four_pts(1,1), four_pts(2,2),temp, 'FontName', 'Times New Roman', 'fontsize', 10)
%     end
    text(mean(four_pts(2:3,1)), four_pts(2,2),temp, 'FontName', 'Times New Roman', 'fontsize', 10)
    
end
scatter(winner(:,1), winner(:,3), 200, '.', 'MarkerEdgeColor',[0.9290, 0.6940, 0.1250])
legend([p1,p2],'Feasibility Boundary','Search Areas', 'fontsize', 12)
legend boxoff
title ('Search Boundaries FW-ST axis', 'fontsize', 16)
xlim([0.9*bounds(1,1) 1.05*bounds(1,2)])
ylim([0.9*bounds(1,5) 1.05*bounds(1,6)])
xlabel('FW [\mum]', 'fontsize', 12)
ylabel('ST [\mum]', 'fontsize', 12)
set(gca, 'FontName', 'Times New Roman')

%% plotting the 3D results
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
title ('Greedy Algorithm Results', 'fontsize', 16)
caxis([min(GreedyData(:,6)) max(GreedyData(:,6))])
colormap(flipud(cmap))
colorbar
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
xlim ([0 1000]*10^-6)
ylim ([0.7 1.5])
title ('Loss as function or pore radius', 'fontsize', 16)
xlabel('r [\mum]')
ylabel('Loss Rate')
set(gca, 'FontName', 'Times New Roman')
legend('Greedy Loss', 'Genetic Loss', 'Analytic Spherical Loss')

%% plotting the loss function against r
subplot ('Position', pos2)
hold on
set(gca, 'FontName', 'Times New Roman')
plot(1:length(GreedyData(:,6)), GreedyData(:,6));
plot(1:length(GeneticData(:,6)), GeneticData(:,6));
title ('Algorithm Results', 'fontsize', 16)
grid on

xlabel('Iteration number', 'fontsize', 12)
ylabel('Loss Rate', 'fontsize', 12)
set(gca, 'FontName', 'Times New Roman')
legend('Greedy Loss', 'Genetic Loss')
