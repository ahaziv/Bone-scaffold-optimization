close all; clear all; clc;
addpath('D:\Documents\2nd degree\Mod Frontier analysis\Curvature Model');

InitData = readtable('SYS.dat','Delimiter','\n','ReadVariableNames',false); 
idx = ismember(InitData.Var1, '-1');
[idxnum,~] = find(idx,2);

opts = delimitedTextImportOptions(...
        'Delimiter', ' ',...
        'ConsecutiveDelimitersRule', 'join','LeadingDelimitersRule','ignore',...
        'MissingRule', 'omitvar','LineEnding','\n','ImportErrorRule','omitvar');

InitData = readtable('SYS.dat',opts); 
% extracting the Elements and Node tables
ElTable = InitData(26:idxnum(1)-1,2:4);
NoTable = InitData(idxnum(1)+7:idxnum(2)-1,12:14);
% converting the tables to arrays
Elements = table2array(ElTable);
Elements = str2double(Elements);
Nodes = table2array(NoTable);
Nodes = str2double(Nodes);


