% this script plots the evolving population sampling
clear all; close all; clc;

FilePath = mfilename('fullpath');
FileName = 'Genetic Data16-12-2019';
idcs = strfind(FilePath,'\');
ParDir = FilePath(1:idcs(end-1)-1);
addpath(strcat(ParDir,'\Data Files'));
RawData = xlsread(FileName,FileName);

Data = RawData((RawData(:,7)==1),1:6);
Data(:,2) = Data(:,2).*Data(:,1);
Data(:,3) = Data(:,3).*Data(:,1);

%% cleaning the data
temp_length = 1;
while temp_length ~= length(Data(:,1))
    Avg_Loss = mean(Data(:,6));
    temp_length = length(Data(:,1));
    Data(Data(:,6)>(50*Avg_Loss),:) = [];
end

%% plotting the 3D results
figure
set(gcf, 'Position', [0,0,960,480]);
cmap = parula;
pos1 = [0.06 0.09 0.43 0.83];
pos2 = [0.55 0.09 0.43 0.83];

c = -log(Data(:,6)-0.9*min(Data(:,6)));
cn = (c-min(c))/(max(c)-min(c));
cn = ceil(cn*size(cmap,1));
cn = max(cn,1);

subplot ('Position', pos2)
set(gca, 'FontName', 'Times New Roman')
scatter3(Data(:,1),Data(:,2),Data(:,3), 150, cmap(cn(:),:), '.');
title ('Genetic Algorithm Results', 'fontsize', 16)
caxis([min(Data(:,6)) max(Data(:,6))])
colormap(flipud(cmap))
colorbar
grid on

xlabel('FW [\mum]', 'fontsize', 12)
ylabel('FD [\mum]', 'fontsize', 12)
zlabel('ST [\mum]', 'fontsize', 12)
set(gca, 'FontName', 'Times New Roman')

%% plotting the 2D results against iterations
subplot ('Position', pos1)
set(gca, 'FontName', 'Times New Roman')
plot(1:length(Data(:,6)), Data(:,6));
title ('Genetic Algorithm Results', 'fontsize', 16)
grid on

xlabel('Iteration number.', 'fontsize', 12)
ylabel('Loss Rate', 'fontsize', 12)
set(gca, 'FontName', 'Times New Roman')