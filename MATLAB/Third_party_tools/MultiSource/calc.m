%% 
clc;
clear;

xray = getGeometry();
xray = calcXrayAngels(xray);
xray = calcVchannelParam(xray);

% 重排成平行束结构后的射线束
