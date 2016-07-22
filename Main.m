% Main ALV-2 Trajectory optimisation Routine
clear all
clc
t1 = cputime;

prompt = {'Launch Altitude (km)','Launch Longitude (deg)','Launch Latitude (deg)', 'Launch Angle (deg)', 'Launch Heading Angle (deg)', 'Target Altitude (km)', 'Second Stage Node Spacing (s)','Third Stage Node Spacing (s)', 'Pre-Pitchover Flight Time (s)', 'Pitchover Angle (deg)','Guess Pitching Angle rad', 'Pitching Angle Optimisation'};
dlg_title = 'Inputs';
num_lines = 1;
defaultans = {'0','153','-27','90','97', '400', '30', '60', '15','1','Auto','NO'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

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
pitchover_angle = pi/180*str2num(answer{10}); % Angle of pitch at the end of the prepitch flight (rad)
Guess = answer{11}; % Either set to 'Auto' for an automatically computed guess, or input a single value (rad) of AoA over the whole 2nd + 3rd stage trajectory
PAO = answer{12}; % Pitching angle optimisation, YES or NO

if strcmp(Guess,'Auto') == 1 % Automatically compute the best guess of constant pitch change rate over range of allowable solutions
  n=1;
  for i = -0.3:0.05:0.200 % guess AoA ranga
    StageDynamics = ALV2Optimiser(icond,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,pitchover_angle,i,'noOpt'); 
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
if strcmp(PAO,'NO') == 1 % No pitching angle optimisation
  [StageDynamics x tspan1 prepitch tspan2 postpitch] = ALV2Optimiser(icond,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,pitchover_angle,Guess,'Opt');
elseif strcmp(PAO,'YES') == 1 % No pitching angle optimisation
  max_v = 0;
  for i = 1:0.5:5
    [temp1 temp2 temp3 temp4 temp5 temp6] = ALV2Optimiser(icond,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,i,Guess,'Opt');
    if temp(end,2) > max_v
      max_v = temp(end,2)
      StageDynamics = temp1;
      x = temp2;
      tspan1 = temp3;
      prepitch = temp4;
      tspan2 = temp5;
      postpitch = temp6;
    end
  end
end

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