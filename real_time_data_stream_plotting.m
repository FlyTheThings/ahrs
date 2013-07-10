%% Real Time Data Stream Plotting Example
function real_time_data_stream_plotting
clc
addpath('quaternion_library');
global AHRS
global Gyroscope Accelerometer Magnetometer
Gyroscope = zeros(1,3)';
Accelerometer = zeros(1,3)';
Magnetometer = zeros(1,3)';
AHRS = MyAHRS();
% AHRS = MadgwickAHRS('SamplePeriod', 1/50, 'Beta', 0.1);
% AHRS = MahonyAHRS('SamplePeriod', 1/100, 'Kp', 0.5);

interfaceObject = udp('0.0.0.0',19,'LocalPort',5555);
interfaceObject.InputBufferSize  = 8192;
%% 
% Setup a figure window and define a callback function for close operation
figureHandle = figure('NumberTitle','off',...
    'Name','Live Data Stream Plot',...
    'Color',[0 0 0],...
    'CloseRequestFcn',{@localCloseFigure,interfaceObject});

%%
% Setup the axes 
axesHandle = axes('Parent',figureHandle,...
    'YGrid','on',...
    'YColor',[0.9725 0.9725 0.9725],...
    'XGrid','on',...
    'XColor',[0.9725 0.9725 0.9725],...
    'Color',[0 0 0]);

xlabel(axesHandle,'Number of Samples');
ylabel(axesHandle,'Value');

%%
% Initialize the plot and hold the settings on
hold on;
plotHandle = plot(axesHandle,0,'-r','LineWidth',1);
plotHandle(2) = plot(axesHandle,0,'-g','LineWidth',1);
plotHandle(3) = plot(axesHandle,0,'-b','LineWidth',1);

%%
% Define a callback function to be executed when desired number of bytes
% are available in the input buffer
interfaceObject.DatagramReceivedFcn = {@localReadAndPlot,plotHandle};

%% 
% Open the interface object
fopen(interfaceObject);
global t
t=0;
global fifo;
fifo = zeros(3,1000);

% Access the 3D World from MATLAB
world=vrworld('my_plane.wrl', 'new');
%world=vrworld('my_box.wrl', 'new');
open(world);
fig=vrfigure(world);
set(fig, 'Viewpoint', 'Far View');
global airpln
airpln=vrnode(world, 'Plane');

%% Implement the bytes available callback
function localReadAndPlot(interfaceObject,x,figureHandle)
global fifo AHRS
%% 
% Read the desired number of data bytes
val=fscanf(interfaceObject,'%f,%d,%f,%f,%f,%d,%f,%f,%f,%d,%f,%f,%f');
if(numel(val)<9)
    disp('Pkt Drop')
    return
end
global Gyroscope Accelerometer Magnetometer airpln
Time = val(1);
Accelerometer = val(3:5)';
Gyroscope = val(7:9)';
if(numel(val)==13)
    Magnetometer = val(11:13)';
end
AHRS.Update(Time, Gyroscope, -Accelerometer, Magnetometer);
if interfaceObject.BytesAvailable == 0
    quaternion = AHRS.Quaternion;
    %euler = quatern2euler(quaternConj(quaternion)) ;

    % Update the plot
    %fifo = [fifo(:,2:end) (180/pi)*euler'];
    %set(figureHandle(1),'Ydata',fifo(1,:));
    %set(figureHandle(2),'Ydata',fifo(2,:));
    %set(figureHandle(3),'Ydata',fifo(3,:));
    %drawnow

    % Rotation setting for the Plane node
    q = quaternConj(quaternion);
    vector = [-q(2)  q(1) q(3)]/ norm(q(1:3));
    theta  = 2 * acos(q(4));

    airpln.rotation=[vector theta];
    % Update the figure
    vrdrawnow
end

%% Implement the close figure callback
function localCloseFigure(figureHandle,x,interfaceObject)

%% 
% Clean up the interface object
fclose(interfaceObject);
delete(interfaceObject);
clear interfaceObject;

%% 
% Close the figure window
delete(figureHandle);
