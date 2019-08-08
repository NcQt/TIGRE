function xray = calcVchannelParam(xray)
% 计算每个角度下探测器通道到虚拟探测器通道的映射
% 得到每个角度下所有探测器通道对应的虚拟探测器通道数
% 计算出正弦图对应的虚拟探测器通道
channel0_xy = [xray.channel_x0(1),xray.channel_y0(1)];
channel1_xy = [xray.channel_x0(2),xray.channel_y0(2)];
origin_xy = xray.origin_xy;
v_distance = xray.VirtualDistance;
detectorVector = xray.detectorVector;
pixelNum1 = xray.channel_num(1);
pixelNum2 = xray.channel_num(2);

vchannelsize = zeros(xray.channel_totalnum,xray.ParallelAngleNum);
vchannelsign = zeros(xray.channel_totalnum,xray.ParallelAngleNum);
vchannelPostion = zeros(xray.channel_totalnum,xray.ParallelAngleNum);
vchannelflag = zeros(xray.channel_totalnum,xray.ParallelAngleNum);
for i = 1:xray.ParallelAngleNum
    parallelangel = xray.ParallelAngle(i);
    v1 = calcVchannel(parallelangel,channel0_xy,origin_xy,v_distance,detectorVector(1,:),pixelNum1);
    v2 = calcVchannel(parallelangel,channel1_xy,origin_xy,v_distance,detectorVector(2,:),pixelNum2);
    vchannelsign(:,i) = [v1.sign;v2.sign];
    vchannelPostion(:,i) = [v1.vchannelPostion;v2.vchannelPostion];
    vchannelsize(:,i) = [v1.channelsize; v2.channelsize];
    %
    for j = 1:length(vchannelPostion(:,i))
        if vchannelsign(j,i) == 1 && vchannelPostion(j,i)>=1 && vchannelPostion(j,i)<=1000
            vchannelflag(j,i) = 1;
        else
            vchannelflag(j,i) = 0;
        end
    end
end

xray.VirtualChannelFlag = vchannelflag;
xray.VirtualChannelPostion = vchannelPostion;
xray.VirtualChannelSize = vchannelsize;

end



function v = calcVchannel(angle,channel0_xy,origin_xy,v_distance,pixelVector,pixelNum)
%% 
% 已知投影角度 v.angle, 探测器通道坐标x,y
% 已知虚拟探测器通道间距
% 求对应虚拟探测器的坐标及通道数
v.angle = angle; % 弧度
v.channel0_xy = channel0_xy;
vpixelsize = 1;

% 计算虚拟探测器通道直线方程
% 虚拟探测器通道中心
% 虚拟探测器通道间隔为1mm

v_beta = v.angle + 90*pi/180;
v_xy = [v_distance*cos(v.angle)+origin_xy(1),v_distance*sin(v.angle)+origin_xy(2)];% 直线方程中心坐标

if abs(v.angle) == 0
    A = [0,1; 1,0];
    B = [v.channel0_xy(2);v_xy(1)];
elseif abs(v.angle) == 0.5*pi
    A = [1,0; 0,1];
    B = [v.channel0_xy(1);v_xy(2)];
else
    A = [-tan(v.angle),1; -tan(v_beta),1];
    B = [v.channel0_xy(2)-v.channel0_xy(1)*tan(v.angle);v_xy(2)-v_xy(1)*tan(v_beta)];
end
% channel0 x y 
v.vchannel0_xy = A\B;% x,y坐标

% 已经知道了坐标，计算通道值
v.distance = sqrt((v.vchannel0_xy(1)-v_xy(1))^2 + (v.vchannel0_xy(2)-v_xy(2))^2) ;
ax = [cos(pi/2+v.angle),sin(pi/2+v.angle)];
bx = [v.vchannel0_xy(1)-v_xy(1),v.vchannel0_xy(2)-v_xy(2)];
v.direction = sign(dot(ax,bx));
v.channel0_num = 500+v.direction*v.distance/vpixelsize ;

% 计算单个探测器通道映射尺寸
v.vchannelSize = dot(pixelVector,[cos(v_beta),sin(v_beta)])/1;
vchannelsign = sign(v.vchannelSize);

v.vchannelPostion = v.channel0_num + (1:1:pixelNum).*v.vchannelSize;
v.vchannelPostion = v.vchannelPostion';

v.sign = ones(size(v.vchannelPostion))*vchannelsign;
v.channelsize = ones(size(v.vchannelPostion))*abs(v.vchannelSize);

v.angle = v.angle/pi*180;
end
