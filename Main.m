% Main ALV-2 Trajectory optimisation Routine
clear all
t1 = cputime;

prompt = {'Launch Altitude (km)','Launch Longitude (deg)','Launch Latitude (deg)', 'Launch Angle (deg)', 'Launch Heading Angle (deg)', 'Target Altitude (km)', 'Second Stage Node Spacing (s)','Third Stage Node Spacing (s)', 'Pre-Pitchover Flight Time (s)', 'Pitchover Angle (deg)','guess Pitching Angle rad'};
dlg_title = 'Inputs';
num_lines = 1;
defaultans = {'0','153','-27','90','97', '400', '30', '60', '15','1','Auto'};
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
SecondStagedt = str2num(answer{7}); 
ThirdStagedt = str2num(answer{8}); 
prepitch_time = str2num(answer{9});
pitchover_angle = pi/180*str2num(answer{10});
Guess = answer{11};

if strcmp(Guess,'Auto') == 1
  % Automatically compute the best guess of constant pitch change rate over range of allowable solutions
  n=1;
  for i = -0.002:0.0002:0.000
    StageDynamics = ALV2Optimiser(icond,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,pitchover_angle,i,'noOpt');
    diff(n,1) = StageDynamics(end,1) - rTarget;
    diff(n,2) = i;
    n = n+1;
  end
  [diff_min,n_min] = min(abs(diff(:,1)));
  Guess = diff(n_min,2);

else
  Guess = str2num(answer{11}); %Guess Pitching Angle, rad
end



% Call Optimisaion Function ----------------------------------------------------

[StageDynamics tspan1 prepitch tspan2 postpitch] = ALV2Optimiser(icond,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,pitchover_angle,Guess,'Opt');


% Polynomial Output ------------------------------------------------------------
% Fits a Polynomial to the optimal trajectory

printf('Polynomial fit of optimal trajectory, coefficients:');
fflush(stdout);
p = polyfit(StageDynamics(:,9),StageDynamics(:,1),5)

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
