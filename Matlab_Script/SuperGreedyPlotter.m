% this script plots the evolving population sampling
clear all; close all; clc;

FlagPlotAnimation = 1;
FlagPlot3D = 1;

FilePath = mfilename('fullpath');
FileName = 'Data23-12-2019 15-48';
idcs = strfind(FilePath,'\');
ParDir = FilePath(1:idcs(end-1)-1);
addpath(strcat(ParDir,'\Data Files'));
RawData = xlsread(FileName);
PauseLength = 0.4;

%% reading the genetic data
Data = RawData((RawData(:,7)==0),1:6);
% GreedyData(:,2) = GreedyData(:,2).*GreedyData(:,1);
% GreedyData(:,3) = GreedyData(:,3).*GreedyData(:,1);
% winner = RawData((RawData(:,7)==1),1:3);
% winner(1,:) = []; 
% winner(:,2) = winner(:,2).*winner(:,1);
% winner(:,3) = winner(:,3).*winner(:,1);

bounds = RawData((RawData(:,7)==2),1:6);
% bounds(:,3) = bounds(:,3).* bounds(:,1);
% bounds(:,4) = bounds(:,4).* bounds(:,2);
% bounds(:,5) = bounds(:,5).* bounds(:,1);
% bounds(:,6) = bounds(:,6).* bounds(:,2);

winner = Data(1,:);
% cleaning the data
temp_length = 1;
while temp_length ~= length(Data(:,1))
    Avg_Loss = mean(Data(:,6));
    temp_length = length(Data(:,1));
    Data(Data(:,6)>(50*Avg_Loss),:) = [];
end

%% plotting FW-FD
pos1 = [0.06 0.09 0.43 0.83];
pos2 = [0.55 0.09 0.43 0.83];
if FlagPlotAnimation
    figure
    set(gcf, 'Position', [0,0,960,480]);
    subplot ('Position', pos1)
    four_pts = [bounds(1,1) bounds(1,3) bounds(1,5)
                bounds(1,1) bounds(1,4) bounds(1,6)
                bounds(1,2) bounds(1,4) bounds(1,6)
                bounds(1,2) bounds(1,3) bounds(1,5)
                bounds(1,1) bounds(1,3) bounds(1,5)];
    plot(four_pts(:,1), four_pts(:,2), 'r','Linewidth',0.75);
    hold on
    title ('Search Boundaries FW-FD axis', 'fontsize', 16)
    xlim([0.9*bounds(1,1) 1.05*bounds(1,2)])
    ylim([0.9*bounds(1,3) 1.05*bounds(1,4)])
    xlabel('FW [\mum]', 'fontsize', 12)
    ylabel('FD/FW [\mum]', 'fontsize', 12)
    set(gca, 'FontName', 'Times New Roman')
    winner_plot1 = scatter(Data(1,1), Data(1,2),300,'.','r');

    subplot ('Position', pos2)
    plot(four_pts(:,1), four_pts(:,3), 'r','Linewidth',0.75);
    hold on
    title ('Search Boundaries FW-ST axis', 'fontsize', 16)
    xlim([0.9*bounds(1,1) 1.05*bounds(1,2)])
    ylim([0.9*bounds(1,5) 1.05*bounds(1,6)])
    xlabel('FW [\mum]', 'fontsize', 12)
    ylabel('ST/FW [\mum]', 'fontsize', 12)
    set(gca, 'FontName', 'Times New Roman')
    winner_plot2 = scatter(Data(1,1), Data(1,3),300,'.','r');

    for ii = 1:length(Data(:,1))
        %% plotting FW - FD/FW
        hold all
        subplot ('Position', pos1)
        scatter(Data(ii,1), Data(ii,2),150,'.','b')
        subplot ('Position', pos2)
        scatter(Data(ii,1), Data(ii,3),150,'.','b')
        if winner(6) > Data(ii,6)
            winner = Data(ii,:);
            delete(winner_plot1)
            delete(winner_plot2)
            subplot ('Position', pos1)
            winner_plot1 = scatter(Data(ii,1), Data(ii,2),300,'.','r');
            subplot ('Position', pos2)
            winner_plot2 = scatter(Data(ii,1), Data(ii,3),300,'.','r');
        end
        four_pts = [bounds(ii,1) bounds(ii,3) bounds(ii,5)
                    bounds(ii,1) bounds(ii,4) bounds(ii,6)
                    bounds(ii,2) bounds(ii,4) bounds(ii,6)
                    bounds(ii,2) bounds(ii,3) bounds(ii,5)
                    bounds(ii,1) bounds(ii,3) bounds(ii,5)];    

        if ii>1
            delete(bound_plot1)
            delete(bound_plot2)
        end
        subplot ('Position', pos1)
        bound_plot1 = plot(four_pts(:,1), four_pts(:,2), 'Color', [0, 0.6, 0],'Linewidth',0.75);
        subplot ('Position', pos2)
        bound_plot2 = plot(four_pts(:,1), four_pts(:,3), 'Color', [0, 0.6, 0],'Linewidth',0.75);
        drawnow
        pause(PauseLength)
        hold off
    end
end
%% plotting the 3D results
if FlagPlot3D
    
    
end


