% Main ALV-2 Trajectory optimisation Routine

% This optimises the rocket trajectory in Earth Centered Rotational coordinates
% ie. to analyse end orbit, look at the velocity and inclination displayed last, rather than heading angle and velocity of rocket, which are Earth relative


clear all
clc
t1 = cputime;

prompt = {'Launch Altitude (km)','Launch Longitude (deg)','Launch Latitude (deg)', 'Launch Angle (deg)', 'Launch Heading Angle (deg)', 'Target Altitude (km)', 'Second Stage Node Spacing (s)','Third Stage Node Spacing (s)', 'Pre-Pitchover Flight Time (s)', 'First-Second Stage Separation (km)','Guess Pitching Angle rad'};
dlg_title = 'Inputs';
num_lines = 1;
defaultans = {'0','153','-27','90','97', '400', '30', '60', '15','42','Auto'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

%coordinates are in rotational geodetic form
icond.h0_prepitch = str2num(answer{1})*1000; % Altitude (m)
icond.xi0_prepitch = pi/180*(str2num(answer{2})); % Longitude (rad)
icond.phi0_prepitch = pi/180*(str2num(answer{3})); % Latitude (rad)
icond.gamma0_prepitch = pi/180*(str2num(answer{4})); % Flight Path Angle (rad)
icond.zeta0_prepitch = pi/180*(str2num(answer{5})); % Heading Angle (rad)
global rTarget
rTarget = str2num(answer{6})*1000; % Target Altitude (m)

% Length of time segments for second + third stage
% Lower lengths mean higher accuracy, but more computation time and potentially convergence issues
SecondStagedt = str2num(answer{7}); % Node separation, second stage (s)
ThirdStagedt = str2num(answer{8}); % Node separation, third stage (s)

prepitch_time = str2num(answer{9}); % Time before pitch is initiated (s)
rTarget_FirstStage = 1000*str2num(answer{10}); % target First-Second Stage Separation (m)
Guess = answer{11}; % Either set to 'Auto' for an automatically computed guess, or input a single value (rad) of AoA over the whole 2nd + 3rd stage trajectory

% Define Vehicle ===============================================================

vehicle.mFirstStage = 480; %(kg) 
vehicle.mFirstStageFuel = 1600;  %(kg) 

vehicle.mSecondStage = 228;
vehicle.mSecondStageFuel = 930;
vehicle.mThirdStage = 40;
vehicle.mThirdStageFuel = 145;
vehicle.mPayload = 18 + 7;

vehicle.mTotal = vehicle.mFirstStage*4 + vehicle.mFirstStageFuel*4 + vehicle.mSecondStage + vehicle.mSecondStageFuel + vehicle.mThirdStage + vehicle.mThirdStageFuel + vehicle.mPayload;

vehicle.mdotFirstStage = 16.39;
vehicle.mdotSecondStage = 3.952;
vehicle.mdotThirdStage = 0.4744;


% Compute Pitchover Angle ======================================================

pitchover_angle = Pitchover(icond,vehicle,rTarget_FirstStage,prepitch_time)

% ==============================================================================


if strcmp(Guess,'Auto') == 1 % Automatically compute the best guess of constant pitch change rate over range of allowable solutions
  n=1;
  for i = -0.3:0.05:0.200 % guess AoA ranga
    StageDynamics = ALV2Optimiser(icond,vehicle,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,pitchover_angle,i,'noOpt'); 
    diff(n,1) = StageDynamics(end,1) - rTarget;
    diff(n,2) = i;
    n = n+1;
  end
  [diff_min,n_min] = min(abs(diff(:,1))); % choose the guess that puts the end of the trajectory closest to target altitude
  Guess = diff(n_min,2);
else
  Guess = str2num(answer{11}); %Input Guess Pitching Angle, rad
end

% Call Optimisaion Function ----------------------------------------------------
[StageDynamics x tspan1 prepitch tspan2 postpitch] = ALV2Optimiser(icond,vehicle,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,pitchover_angle,Guess,'Opt');

Omega_E = 7.2921e-5 ; % rotation rate of the Earth rad/s

v = StageDynamics(end,2);

zeta = StageDynamics(end,8);

rEarth = 6.3674447e6;  %(m) radius of earth

r = StageDynamics(end,1) + rEarth ;

phi = StageDynamics(end,7);

fprintf('Orbital Velocity and inclination');

vexo = sqrt((v*sin(zeta))^2 + (v*cos(zeta) + r*Omega_E*cos(phi))^2) %Change coordinate system when exoatmospheric, add velocity component from rotation of the Earth

inc = asin(v*sin(pi-zeta)/vexo)  % orbit inclination angle



% Guidance Output ------------------------------------------------------------
% Numerical Differentiation of pitching angle, to obtain first stage pitch rate
gammadot1(:,1) = tspan2(1:end-1);
for i = 1:length(tspan2)-1
  gammadot1(i,2) = (postpitch(i+1,4) - postpitch(i,4))/(tspan2(i+1)-tspan2(i));
end

% Fits a Polynomial to the pitching angle time history
printf('Polynomial fit of optimal pitching angle time history, coefficients:');
p = polyfit(StageDynamics(:,9),(StageDynamics(:,4)+StageDynamics(:,4)),5)
fflush(stdout);


% Calculate Runtime
t2 = cputime;
runtime = t2-t1

% Plotting ---------------------------------------------------------------------
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



function par_vend = parfun


end