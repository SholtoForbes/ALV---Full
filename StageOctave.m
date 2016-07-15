% Created by Sholto Forbes 15/7/16
% Single shooting routine for the calculation of ALV-2 orbital trajectories

clc; clear;
clear all;
t1 = cputime;

prompt = {'Launch Altitude (km)','Launch Longitude (deg)','Launch Latitude (deg)', 'Launch Angle (deg)', 'Launch Heading Angle (deg)', 'Target Altitude (km)', 'Second Stage Segment length (s)','Third Stage Segment length (s)'};
dlg_title = 'Inputs';
num_lines = 1;
defaultans = {'0','153','-27','90','97', '400', '30', '60'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

h0_prepitch = str2num(answer{1})*1000; % Altitude (m)
xi0_prepitch = pi/180*(str2num(answer{2})); % Longitude (rad)
phi0_prepitch = pi/180*(str2num(answer{3})); % Latitude (rad)
gamma0_prepitch = pi/180*(str2num(answer{4})); % Flight Path Angle (rad)
zeta0_prepitch = pi/180*(str2num(answer{5})); % Heading Angle (rad)
global rTarget
rTarget = str2num(answer{6})*1000; % Target Altitude (m)

% Length of time segments for second + third stage
% Lower lengths mean higher accuracy, but more computation time and potentially convergence issues
SecondStagedt = str2num(answer{7}); 
ThirdStagedt = str2num(answer{8}); 

   
mFirstStage = 480; %(kg)  %Total lift-off mass
mFirstStageFuel = 1600;  %(kg)  %mass of the fuel
global mSecondStage
mSecondStage = 228;
mSecondStageFuel = 930;
mThirdStage = 40;
mThirdStageFuel = 145;
mPayload = 18 + 7;

mTotal = mFirstStage*4 + mFirstStageFuel*4 + mSecondStage + mSecondStageFuel + mThirdStage + mThirdStageFuel + mPayload;

mdotFirstStage = 16.39;
mdotSecondStage = 3.952;
mdotThirdStage = 0.4744;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

v0_prepitch = 0;  %Rocket starts stationary
m0_prepitch = mTotal;  %Rocket starts full of fuel


phase = 'prepitch';
tspan1 = [0:15]; % prepitch time 
prepitch0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, xi0_prepitch, phi0_prepitch, zeta0_prepitch];

[prepitch] = lsode(@(prepitch,t) rocketDynamics(prepitch,0,phase), prepitch0, tspan1);  
 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Post-Pitchover Simulation                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

phase = 'postpitch';
tspan2 = [15:mFirstStageFuel*4/(mdotFirstStage*4)]; % run until fuel ends
tspan2(end) = mFirstStageFuel*4/(mdotFirstStage*4);

postpitch0 = [prepitch(end,1), prepitch(end,2), prepitch(end,3), 1.55334, 0, prepitch(end,6), prepitch(end,7), prepitch(end,8)];
dalphadt = 0; % gravity turn, AoA=0

[postpitch] = lsode(@(postpitch,t) rocketDynamics(postpitch,dalphadt,phase), postpitch0, tspan2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Stage 2 + 3                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

global initialdynamics
initialdynamics = [postpitch(end,1) postpitch(end,2) postpitch(end,3)-mFirstStage*4 postpitch(end,4) postpitch(end,5) postpitch(end,6) postpitch(end,7) postpitch(end,8)];

% Define Control Time Spans ----------------------------------------------------
%% calculate the second stage runtime, break into segments of constant control
global secondstagetimenodes
secondstagetimenodes = [tspan2(end):SecondStagedt:tspan2(end) + mSecondStageFuel/mdotSecondStage];
secondstagetimenodes(end) = tspan2(end) + mSecondStageFuel/mdotSecondStage;
global n2
n2 = floor(mSecondStageFuel/mdotSecondStage/SecondStagedt);

%% calculate the third stage runtime, break into segments of constant control
global thirdstagetimenodes
thirdstagetimenodes = [secondstagetimenodes(end):ThirdStagedt:secondstagetimenodes(end) + mThirdStageFuel/mdotThirdStage];
thirdstagetimenodes(end) = secondstagetimenodes(end) + mThirdStageFuel/mdotThirdStage;
n3 = floor(mThirdStageFuel/mdotThirdStage/ThirdStagedt);

% Main Second + Third Stage Function -------------------------------------------
%%%%% takes input of dAlpha/dt vector, outputs system dynamics

function [StageDynamics] = Stages(dalphadt)

%Begin Second Stage ------------------------------------------------------------
phase = 'secondstage';

global initialdynamics
global secondstagetimenodes

dynamics0 = initialdynamics; % initialize dynamics

StageDynamics = [dynamics0 secondstagetimenodes(1)]; % add time storage

% Break down each segment of constant dAlpha/dt into 1 second intervals
for i = 1:length(secondstagetimenodes)-1 

tspan = secondstagetimenodes(i):secondstagetimenodes(i+1);

% Calculate dynamics over each segment
[dynamics] = lsode(@(dynamics,t) rocketDynamics(dynamics,dalphadt(i),phase), dynamics0, tspan);

dynamics0 = dynamics(end,:);

% Store dynamics
StageDynamics = [StageDynamics;[dynamics transpose(tspan)]];

end

% no second stage nodes
global nSecondStage
nSecondStage = length(StageDynamics(:,1));

% Stage 2-> 3 mass change
global mSecondStage
dynamics0 = [dynamics(end,1) dynamics(end,2) dynamics(end,3)-mSecondStage dynamics(end,4) dynamics(end,5) dynamics(end,6) dynamics(end,7) dynamics(end,8)];

%Begin Third Stage -------------------------------------------------------------
phase = 'thirdstage';

global thirdstagetimenodes

% Break down each segment of constant dAlpha/dt into 1 second intervals
for i = 1:length(thirdstagetimenodes)-1

tspan = thirdstagetimenodes(i):thirdstagetimenodes(i+1);

% Calculate dynamics over each segment
[dynamics] = lsode(@(dynamics,t) rocketDynamics(dynamics,dalphadt(i+length(secondstagetimenodes)-1),phase), dynamics0, tspan);

dynamics0 = dynamics(end,:); 

% Store dynamics
StageDynamics = [StageDynamics;[dynamics transpose(tspan)]];
end

endfunction

% Functions which call main dynamics function ----------------------------------
%% these simply return specific outputs that the SQP function requires

% Function to call velocity, used to maximise 
function vend = velocityend(dalphadt)

[StageDynamics] = Stages(dalphadt);

vend = -StageDynamics(end,2);

endfunction

% Function to call altitude - target altitude, to constrain end altitude
function constraints = const(dalphadt)

[StageDynamics] = Stages(dalphadt);

hend = StageDynamics(end,1);

global rTarget
constraints = [hend-rTarget];

endfunction

% function to define maximum altitude, this helps to guide the rocket into the correct orbital position
function ineq = inequalities(dalphadt)

[StageDynamics] = Stages(dalphadt);
global rTarget
ineq = -StageDynamics(:,1)+rTarget;

endfunction

% GUESS 
%This is important, must be same no. nodes as required output, must be near the expected answer to within reason

x0 = zeros(1,n2+n3)-0.001;

  
% SQP --------------------------------------------------------------------------
% This is the main optimisation routine utilising Sequential Quadratic Programming
[x, obj, inform, iter, nf, lambda] = sqp (x0, @velocityend, @const, @inequalities, -0.002, 0.002);

%Forward Simulation
[StageDynamics] = Stages(x);


t2 = cputime;
runtime = t2-t1

% Plotting
global nSecondStage
figure(1)

subplot(5,1,1)
hold on
plot(tspan1,prepitch(:,1)/1000,'LineWidth',1.5,'Color','b')
plot(tspan2,postpitch(:,1)/1000,'LineWidth',1.5,'Color','r')
plot(StageDynamics(1:nSecondStage,9),StageDynamics(1:nSecondStage,1)/1000,'LineWidth',1.5,'Color','g')
plot(StageDynamics(nSecondStage:end,9),StageDynamics(nSecondStage:end,1)/1000,'LineWidth',1.5,'Color','y')
ylabel('Altitude (km)');
subplot(5,1,2)
hold on
plot(tspan1,180/pi*prepitch(:,4),'LineWidth',1.5,'Color','b')
plot(tspan2,180/pi*postpitch(:,4),'LineWidth',1.5,'Color','r')
plot(StageDynamics(1:nSecondStage,9),180/pi*StageDynamics(1:nSecondStage,4),'LineWidth',1.5,'Color','g')
plot(StageDynamics(nSecondStage:end,9),180/pi*StageDynamics(nSecondStage:end,4),'LineWidth',1.5,'Color','y')
ylabel('Trajectory Angle (deg)');
subplot(5,1,3)
hold on
plot(tspan1,prepitch(:,2),'LineWidth',1.5,'Color','b')
plot(tspan2,postpitch(:,2),'LineWidth',1.5,'Color','r')
plot(StageDynamics(1:nSecondStage,9),StageDynamics(1:nSecondStage,2),'LineWidth',1.5,'Color','g')
plot(StageDynamics(nSecondStage:end,9),StageDynamics(nSecondStage:end,2),'LineWidth',1.5,'Color','y')
ylabel('Velocity (m/s)');
subplot(5,1,4)
hold on
plot(tspan1,prepitch(:,3),'LineWidth',1.5,'Color','b')
plot(tspan2,postpitch(:,3),'LineWidth',1.5,'Color','r')
plot(StageDynamics(1:nSecondStage,9),StageDynamics(1:nSecondStage,3),'LineWidth',1.5,'Color','g')
plot(StageDynamics(nSecondStage:end,9),StageDynamics(nSecondStage:end,3),'LineWidth',1.5,'Color','y')
ylabel('Mass (kg)');
subplot(5,1,5)
hold on
plot(tspan1,180/pi*prepitch(:,5),'LineWidth',1.5,'Color','b')
plot(tspan2,180/pi*postpitch(:,5),'LineWidth',1.5,'Color','r')
plot(StageDynamics(1:nSecondStage,9),180/pi*StageDynamics(1:nSecondStage,5),'LineWidth',1.5,'Color','g')
plot(StageDynamics(nSecondStage:end,9),180/pi*StageDynamics(nSecondStage:end,5),'LineWidth',1.5,'Color','y')
ylabel('Angle of Attack (deg)');
xlabel('Time (s)');

