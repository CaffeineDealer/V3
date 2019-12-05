%% Written by Yavar Korkian - sometimes 2019
% Last modified on 05/12/2019
% This script uses custom written functions including: polarplot, make_row, and rayleighz2
clear all; clear cache
clc
cd E:\MT_MST\mu_wipped\results
type = 'Translation'; % Translation, Spiral, Deformation
load(sprintf('%s.mat',type))
%% Vector Direction
v3Vecsu = horzcat(su.vecV3);
v3Vecmu = horzcat(mu.vecV3);
figure
for i = 1:size(v3Vecsu,2)
    [thV3(i,1),pV3(i,1),rV3(i,1)] = polarplot(v3Vecsu(:,i));
    [thV3(i,2),pV3(i,2),rV3(i,2)] = polarplot(v3Vecmu(:,i));
end
thetaV3 = rad2deg(circ_mean(thV3));

mtVecsu = horzcat(su.vecMT);
mtVecmu = horzcat(mu.vecMT);
for i = 1:size(mtVecsu,2)
    [thMT(i,1),pMT(i,1),rMT(i,1)] = polarplot(mtVecsu(:,i));
    [thMT(i,2),pMT(i,2),rMT(i,2)] = polarplot(mtVecmu(:,i));
end
thetaMT = rad2deg(circ_mean(thMT));
%% Correlation Coefficients
[rv3PD pv3PD] = corrcoef(thV3(:,1),thV3(:,2));
[rmtPD pmtPD] = corrcoef(thMT(:,1),thMT(:,2));
%% Plots
figure;
scatter(thV3(:,1),thV3(:,2),'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
scatter(thMT(:,1),thMT(:,2),'MarkerEdgeColor','g','MarkerFaceColor','g')
refline(1,0);
title(sprintf('%s| PD (V3 = %.2f*, MT = %.2f*)',type,rv3PD(1,2),rmtPD(1,2)));
xlabel 'SU (rad)'; ylabel 'MU (rad)';
legend({'V3','MT'})
