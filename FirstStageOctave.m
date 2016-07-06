

clc; clear;
clear all;
addpath TrajOpt-master

   
mRocket = 27000; %(kg)  %Total lift-off mass
mFuel = 0.8*mRocket;  %(kg)  %mass of the fuel
mSpartan = 8755.1;
mTotal = mSpartan + mRocket;
mEmpty = mRocket-mFuel;  %(kg)  %mass of the rocket (without fuel)
global Tmax
Tmax = 460000;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
h0_prepitch = 0;  %Rocket starts on the ground
v0_prepitch = 0;  %Rocket starts stationary
m0_prepitch = mTotal;  %Rocket starts full of fuel
gamma0_prepitch = 1.5708;

phase = 'prepitch';
tspan = [1:15];
y0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch];

[y] = lsode(@(y,t) rocketDynamics(y,0,phase), y0, tspan);  



phase = 'postpitch';
Tratio = .94;
tspan = [0:(y(end,3)-(mEmpty+mSpartan))/(Tratio*60*Tmax/200000)];
 postpitch0 = [y(end,1), y(end,2), y(end,3), 1.55334];
[postpitch] = lsode(@(postpitch,t) rocketDynamics(postpitch,Tratio*Tmax,phase), postpitch0, tspan);

y
postpitch
postpitch(end,4)






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
