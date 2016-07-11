

clc; clear;
clear all;

   
   
mFirstStage = 480; %(kg)  %Total lift-off mass
mFirstStageFuel = 1600;  %(kg)  %mass of the fuel
mSecondStage = 228;
mSecondStageFuel = 930;
mThirdStage = 40;
mThirdStageFuel = 145;
mPayload = 18 + 7;

mTotal = mFirstStage*4 + mFirstStageFuel*4 + mSecondStage + mSecondStageFuel + mThirdStage + mThirdStageFuel + mPayload;

mdotFirstStage = 16.39;
mdotSecondStage = 3.952;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = 0;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
m0_prepitch = mTotal;  %Rocket starts full of fuel
gamma0_prepitch = 1.5708;

phase = 'prepitch';
tspan = [0:15];
y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0];

[y] = lsode(@(y,t) rocketDynamics(y,0,phase), y0, tspan);  
 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Post-Pitchover Simulation                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

phase = 'postpitch';
tspan = [15:mFirstStageFuel*4/(mdotFirstStage*4)];
postpitch0 = [y(end,1), y(end,2), y(end,3), 1.55334, 0];
dalphadt = 0;
[postpitch] = lsode(@(postpitch,t) rocketDynamics(postpitch,dalphadt,phase), postpitch0, tspan);

y;
postpitch;
postpitch(end,4);
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Stage 2                                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

global initialdynamics
initialdynamics = [postpitch(end,1) postpitch(end,2) postpitch(end,3)-mFirstStage*4 postpitch(end,4) postpitch(end,5)];


global timenodes
timenodes = [0:30:mSecondStageFuel/mdotSecondStage];
timenodes(end) = mSecondStageFuel/mdotSecondStage
n = floor(mSecondStageFuel/mdotSecondStage/30);



function [v23 alphaend SecondStageDynamics] = SecondStage(dalphadt)

phase = 'secondstage';

global initialdynamics
global timenodes

dynamics0 = initialdynamics;

SecondStageDynamics(1,:) = dynamics0;

for i = 1:length(timenodes)-1

tspan = timenodes(i):timenodes(i+1);

[dynamics] = lsode(@(dynamics,t) rocketDynamics(dynamics,dalphadt(i),phase), dynamics0, tspan);

dynamics0 = dynamics(end,:);

SecondStageDynamics(i+1,:) = dynamics(end,:);

end

v23 = -SecondStageDynamics(end,2); % negative, so that minimiser gives maximum

alphaend = SecondStageDynamics(end,5);

endfunction

function v23 = velocity23(dalphadt)

[v23 alphaend SecondStageDynamics] = SecondStage(dalphadt);

endfunction

function alphaend = alpha23(dalphadt)

[v23 alphaend SecondStageDynamics] = SecondStage(dalphadt);

endfunction


x0 = zeros(1,n);
[x, obj, info, iter, nf, lambda] = sqp (x0, @velocity23, @alpha23, [], -0.01, 0.01);

[v23 alphaend SecondStageDynamics] = SecondStage(x);





%
%% Forward Simulation ======================================================
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%                        Pre-Pitchover Simulation                         %
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
%% Note this doesnt work for AoA control
%
%f_h0_prepitch = 0;  %Rocket starts on the ground
%f_v0_prepitch = 0;  %Rocket starts stationary
%f_m0_prepitch = mTotal;  %Rocket starts full of fuel
%f_gamma0_prepitch = deg2rad(90);
%
%phase = 'prepitch';
%f_tspan = [0 15];
%f_y0 = [f_h0_prepitch, f_v0_prepitch, f_m0_prepitch, f_gamma0_prepitch];
%% [f_t_prepitch, f_y_prepitch] = ode45(@(f_t,f_y) rocketDynamics(f_y,Tmax,phase), f_tspan, f_y0);
%[f_t_prepitch, f_y_prepitch] = ode45(@(f_t,f_y) rocketDynamics(f_y,0,phase), f_tspan, f_y0);
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%                        Post-Pitchover Simulation                         %
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
%f_h0 = f_y_prepitch(end,1);  %Rocket starts on the ground
%f_v0 = f_y_prepitch(end,2);  %Rocket starts stationary
%f_m0 = f_y_prepitch(end,3);  %Rocket starts full of fuel
%f_gamma0 = deg2rad(89);    % pitchover 
%
%phase = 'postpitch';
%f_tspan = [0 t(end)];
%f_y0 = [f_h0, f_v0, f_m0, f_gamma0];
%[f_t, f_y] = ode45(@(f_t,f_y) rocketDynamics(f_y,ControlFunction(f_t,t,u),phase), f_tspan, f_y0);
%
%
%
%
%% Plotting
%
%figure(120);
%subplot(2,3,1);
%hold on
%plot(t,x(1,:)/1000)
%plot(f_t,f_y(:,1)/1000)
%xlabel('time (s)')
%ylabel('height (km)')
%subplot(2,3,2);
%plot(t,x(3,:))
%xlabel('time (s)')
%ylabel('mass (kg)')
%subplot(2,3,3);
%plot(t,x(2,:))
%xlabel('time (s)')
%ylabel('velocity (m/s)')
%subplot(2,3,4);
%plot(t,x(4,:))
%xlabel('time (s)')
%ylabel('trajectory angle (rad)')
%subplot(2,3,5);
%plot(t,u)
%xlabel('time (s)')
%ylabel('Control')
