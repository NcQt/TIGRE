function VirtualParam = xrayRebin(xray)
%% 
% 输入角度，几何参数，探测器通道
% 输出对应的焦点位置，相邻焦点

% 1.计算该通道对应的角度范围
% 2.查找焦点编号
% 3.直线相交计算实际焦点位置

timeID = tic;

 xray.ParallelAngle = -135:1:45;
 xray.ParallelAngleNum = length(xray.ParallelAngle);
% 虚拟探测器参数
% parallel.VirtualPixelSize = 1;  % 虚拟探测器通道间距为2mm；
% parallel.VirtualPixelNum = 1000; % 虚拟探测器通道数为500；

virtualNum = 1;
for k = 1:xray.ParallelAngleNum
    
    for i = 1:xray.channel_totalnum
        % 1.计算该通道对应的角度范围
        
        if xray.ParallelAngle(k)>= xray.ChannelThetaRange(i,1)...
                && xray.ParallelAngle(k)<= xray.ChannelThetaRange(i,2);
            % 2.查找焦点编号
            VirtualParam(virtualNum).channelNum = i;
            VirtualParam(virtualNum).angelNum = k;
            for j = 1:xray.source_totalnum-1
                %disp([num2str(xray.ParallelAngle(k)),' ',num2str(xray.theta(j,i)),' ',num2str(xray.theta(j+1,i))]);
                if xray.ParallelAngle(k)<= xray.theta(j,i)...
                        && xray.ParallelAngle(k)>= xray.theta(j+1,i);
                    theta1 = abs(xray.theta(j,i)-xray.ParallelAngle(k));
                    theta2 = abs(xray.theta(j+1,i)-xray.ParallelAngle(k));
                    w = theta1/(theta1+theta2);
                    sourceParam.sourceNum = [j,j+1];
                    sourceParam.sourceWedge = [1-w,w];
                    VirtualParam(virtualNum).sourceParam = sourceParam;
                    break;
                end
            end
%             3.计算虚拟通道位置及对应权重
%             每个角度下计算每排探测器起始通道的虚拟探测器通道数，从而标定出所有的通道对应的虚拟通道
            
            ChannelPostion = xray.VirtualChannelPostion(i,k);
            ChannelFlag = xray.VirtualChannelFlag(i,k);
            ChannelSize = xray.VirtualChannelSize(i,k);
%             虚拟探测器通道编号

            chParam = getChannelPostionNum(ChannelPostion,ChannelSize);
            VirtualParam(virtualNum).chParam = chParam;
            virtualNum = virtualNum + 1;
        end 
    end
end
elapsedTime = toc(timeID);
disp(['The elapsedTime for xrayRebin is ' num2str(elapsedTime)]);
end


function chParam = getChannelPostionNum(ChannelPostion,ChannelSize)
% channel postion num
ChannelPostionUp = ChannelPostion + ChannelSize;
ChannelPostionLow = ChannelPostion;
vChannelPostionNum = floor(ChannelPostionLow):1:floor(ChannelPostionUp);

% channel wedge
if length(vChannelPostionNum)>1
    vChannelPostion1 = [vChannelPostionNum(2:end),ChannelPostionUp];
    vChannelPostion2 = [ChannelPostionLow,vChannelPostionNum(2:end)];
else
    vChannelPostion1 = ChannelPostionUp;
    vChannelPostion2 = ChannelPostionLow; 
end
vChannelWegde = vChannelPostion1 - vChannelPostion2;

chParam.ChannelPostionNum = vChannelPostionNum;
chParam.ChannelPostionWedge = vChannelWegde;
end