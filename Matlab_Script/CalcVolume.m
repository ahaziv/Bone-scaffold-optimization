function [TotVolume] = CalcVolume(FW, FD, ST)
% this function calculates the void space in the scaffold model by
% substituting the scaffold volume from the total unit cube volume.
% the scaffold volume is calculated as the volume of 2 cylinders with a radius
% of FW/2 and hight of FD. from this the function substitutes the
% intersection volume wihch is calculated by a volumetric integration.
% FW = 0.1*10^-3; FD=0.2*10^-3; ST = 0.08*10^-3; just to check stuff
syms x
R = FW/2;
h = ST/2;

CubeVol = FD*FD*2*ST;
FilamentVol = 2*FD*pi()*(FW/2)^2;
% fun = @(x) pi-2*asin((h-sqrt(R^2-x.^2))/R+1)-sin(2*asin((h-sqrt(R^2-x.^2))/R+1));
% IntersecVol = R^2*integral(fun,0,sqrt(R^2-h^2));

fun = @(x) 4*sqrt(R^2-(x+h).^2).*(R^2-(x-h).^2);
IntersecVol = integral(fun,-R+h,R-h);
materialVol = (FilamentVol-IntersecVol)*10^9;
TotVolume = CubeVol - FilamentVol + IntersecVol;