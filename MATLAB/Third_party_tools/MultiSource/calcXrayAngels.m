function xray = calcXrayAngels(xray)
%% 1. 计算焦点对应的探测器通道的所有角度
xray.theta = zeros(xray.source_totalnum,xray.channel_totalnum);
xray.length = zeros(xray.source_totalnum,xray.channel_totalnum);
for i = 1:xray.source_totalnum
    for j = 1:xray.channel_totalnum
        x1 = xray.channel_xy(j,1) - xray.source_xy(i,1);
        y1 = xray.channel_xy(j,2) - xray.source_xy(i,2);
        x2 = 1;
        y2 = 0;
        xray.length(i,j) = sqrt(x1*x1+y1*y1); 
        xray.theta(i,j) = abs(acos((x1*x2+y1*y2)/(xray.length(i,j)*sqrt(x2*x2+y2*y2))))*180/pi;
        if y1<0
           xray.theta(i,j) = -1*xray.theta(i,j);
        end
    end
end

% theta range
a = min(xray.theta);
a = a';
b = max(xray.theta);
a(:,2) = b';
xray.ChannelThetaRange = a;

a1 = min(xray.theta,[],2);
b1 = max(xray.theta,[],2);
a1(:,2) = b1;
xray.SourceThetaRange = a1;

% 画出射线角度
f = figure;
f.Position = [595 100 560 540];
plot(xray.theta');
xlim([0,xray.channel_totalnum]);
xlabel('channel num');
ylabel('xray angle / °');

end
