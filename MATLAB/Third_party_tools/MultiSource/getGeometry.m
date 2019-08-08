function xray = getGeometry()
%% geometry 
% 2d 
% detector channels 1mm*1024channels,rectangular detector struct
% X-ray tube 128 focac spots ,2 tubes ,rectangular geometry

% 重建中心为原点，以下坐标系长度单位均为mm
xray.origin_xy = [0 0];
xray.reconsize = 500 ;% 单位mm， 重建范围为500mm半径的圆

xray.tube_num = 2;
xray.source_x0 = [-500 -496]; % 每根射线管的初始焦点x坐标
xray.source_y0 = [-500 500]; % 每根射线管的初始焦点y坐标
xray.source_delt = [4 4]; % 每根射线管的焦点间距4mm
xray.source_num = [250 249]; % 每根射线管工作的焦点个数
xray.source_totalnum = sum(xray.source_num);
xray.tubeVector = [0,1;1,0]; % 射线管的方向向量

xray.detector_num = 2;
xray.channel_x0 = [-499 500]; % 每段探测器初始通道的x坐标
xray.channel_y0 = [-500 -500]; % 每段探测器初始通道的y坐标
xray.channel_delt = [1 1]; % 探测器通道的间距1mm
xray.channel_num = [999 1000]; % 每段探测器的通道数
xray.channel_totalnum = sum(xray.channel_num);
xray.detectorVector = [1,0;0,1]; % 探测器的方向向量

% 平行束参数
xray.ParallelAngle = (-135:1:45)./180*pi;
xray.ParallelAngleNum = length(xray.ParallelAngle);
 
% 虚拟探测器参数
xray.VirtualPixelSize = 1;  % 虚拟探测器通道间距为1mm；
xray.VirtualPixelNum = 1000; % 虚拟探测器通道数为1000；
xray.VirtualDistance = 510; % 虚拟探测器到重建中心的距离为510mm;
xray.VirtualAngle = xray.ParallelAngle + 0.5*pi; % 弧度制，虚拟探测器的角度
xray.VirtualVector = [cos(xray.VirtualAngle),sin(xray.VirtualAngle)];



%% 射线管焦点坐标
xray.source_xy = zeros([sum(xray.source_num),2]);
for i = 1:xray.tube_num
    if i == 1
        for j = 1:xray.source_num(i)
            xray.source_xy(j,1) = xray.source_x0(i);
            xray.source_xy(j,2) = xray.source_y0(i)+(j-1)*xray.source_delt(i);
        end
    else
        for j = 1:xray.source_num(i)
            xray.source_xy(j+xray.source_num(1),1) = xray.source_x0(i)+(j-1)*xray.source_delt(i);
            xray.source_xy(j+xray.source_num(1),2) = xray.source_y0(i);
        end
    end
end

% 探测器通道坐标
xray.channel_xy = zeros([sum(xray.channel_num),2]);
for i = 1:xray.detector_num
    if i == 1
        for j = 1:xray.channel_num(i)
            xray.channel_xy(j,1) = xray.channel_x0(i)+(j-1)*xray.channel_delt(i);
            xray.channel_xy(j,2) = xray.channel_y0(i);
        end
    else
        for j = 1:xray.channel_num(i)
            xray.channel_xy(j+xray.channel_num(1),1) = xray.channel_x0(i);
            xray.channel_xy(j+xray.channel_num(1),2) = xray.channel_y0(i)+(j-1)*xray.channel_delt(i);
        end
    end
end

%% plot
plotflag = 1;
savegif = 1;

if plotflag
    if savegif
        savename = 'test.gif';
        showGeometry(xray,savegif,savename);
    else
        showGeometry(xray);
    end
end

end


function showGeometry(varargin)
%%
if nargin < 3
    xray = varargin{1};
    savegif = false;
else
    xray = varargin{1};
    savegif = varargin{2};
    savename = varargin{3};
end

%% plot
f = figure;
f.Position = [35 100 560 540];
if savegif == true
    filename = savename;
end

for i = 1:6:xray.source_totalnum
    plot(xray.channel_xy(:,1),xray.channel_xy(:,2));
    hold on;
    plot(xray.source_xy(:,1),xray.source_xy(:,2),'-*');
    rectangle('Position',[-500 -500 1000 1000]);
    rectangle('Position',[-500 -500 1000 1000],'Curvature',[1 1]);
    for j = 1:30:xray.channel_totalnum
        x1 = xray.source_xy(i,1);
        y1 = xray.source_xy(i,2);
        x2 = xray.channel_xy(j,1);
        y2 = xray.channel_xy(j,2);
        plot([x1 x2],[y1 y2],'r');
    end
    hold off;
    drawnow;
    if savegif
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1;
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    end  
end
end