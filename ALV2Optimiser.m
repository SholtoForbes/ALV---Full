% Created by Sholto Forbes 15/7/16
% Single shooting routine for the calculation of ALV-2 orbital trajectories

function [StageDynamics x tspan1 prepitch tspan2 postpitch] = ALV2Optimiser(icond,rTarget,SecondStagedt,ThirdStagedt,prepitch_time,pitchover_angle,Guess, optim)

h0_prepitch = icond.h0_prepitch;
xi0_prepitch = icond.xi0_prepitch;
phi0_prepitch = icond.phi0_prepitch;
gamma0_prepitch = icond.gamma0_prepitch;
zeta0_prepitch = icond.zeta0_prepitch;
   
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
printf('Pre-Pitch Running    ')
fflush(stdout);

v0_prepitch = 0;  %Rocket starts stationary
m0_prepitch = mTotal;  %Rocket starts full of fuel

phase = 'prepitch';
tspan1 = [0:prepitch_time]; % prepitch time 
prepitch0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, xi0_prepitch, phi0_prepitch, zeta0_prepitch];

[prepitch] = lsode(@(prepitch,t) rocketDynamics(prepitch,t,0,phase,tspan1), prepitch0, tspan1);  
 

printf('Complete    ')
fflush(stdout);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Post-Pitchover Simulation                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
printf('Post-Pitch Running    ')
fflush(stdout);

phase = 'postpitch';
tspan2 = [15:mFirstStageFuel*4/(mdotFirstStage*4)]; % run until fuel ends
tspan2(end) = mFirstStageFuel*4/(mdotFirstStage*4);

postpitch0 = [prepitch(end,1), prepitch(end,2), prepitch(end,3), pi/2-pitchover_angle, 0, prepitch(end,6), prepitch(end,7), prepitch(end,8)];
dalphadt = 0; % gravity turn, AoA=0

[postpitch] = lsode(@(postpitch,t) rocketDynamics(postpitch,t,dalphadt,phase,tspan2), postpitch0, tspan2);

printf('Complete    ')
fflush(stdout);
 
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

% GUESS 
%This is important, must be same no. nodes as required output, must be near the expected answer to within reason

x0 = zeros(1,n2+n3)+Guess;

  
switch optim % allows routine to be run without optimisation for initial guess checking
  case 'Opt'
     printf('Single Shooting Optimiser Running    ')
     fflush(stdout);
     
    % SQP --------------------------------------------------------------------------
    % This is the main optimisation routine utilising Sequential Quadratic Programming
    [x, obj, inform, iter, nf, lambda] = sqp (x0, @cost, @const, @inequalities, -0.4, 0.4, 100, 1e-04);
    
    %Forward Simulation
    [StageDynamics] = Stages(x);

  case 'noOpt'
    printf('Second and Third Stage Running    ')
    fflush(stdout);
    x = x0;
    [StageDynamics] = Stages(x);
  end
  printf('Complete    ')
  fflush(stdout);

endfunction

function [StageDynamics] = Stages(dalphadt)

%Begin Second Stage ------------------------------------------------------------
phase = 'secondstage';

global initialdynamics
global secondstagetimenodes

dynamics0 = initialdynamics; % initialize dynamics

StageDynamics = [dynamics0 secondstagetimenodes(1)]; % add time storage

% Break down each segment of constant dAlpha/dt into 1 second intervals
for i = 1:length(secondstagetimenodes)-1 

  tspan = [secondstagetimenodes(i) secondstagetimenodes(i+1)];

  % Calculate dynamics over each segment
  %lsode_options ("minimum step size", .001);
  [dynamics] = lsode(@(dynamics,t) rocketDynamics(dynamics,t,dalphadt(i),phase,tspan), dynamics0, tspan);

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

  tspan = [thirdstagetimenodes(i) thirdstagetimenodes(i+1)];
  % Calculate dynamics over each segment
  [dynamics] = lsode(@(dynamics,t) rocketDynamics(dynamics,t,dalphadt(i+length(secondstagetimenodes)-1),phase,tspan), dynamics0, tspan);

  dynamics0 = dynamics(end,:); 

  % Store dynamics
  StageDynamics = [StageDynamics;[dynamics transpose(tspan)]];

end

printf('End Altitude Achieved')
disp(StageDynamics(end,1));
fflush(stdout);
endfunction

% Functions which call main dynamics function ----------------------------------
%% these simply return specific outputs that the SQP function requires

% Function to call velocity, used to maximise 
function vend = cost(dalphadt)
  [StageDynamics] = Stages(dalphadt);
  hend = StageDynamics(end,1);
  global rTarget
  vend = -StageDynamics(end,2) ;
endfunction

% Function to call altitude - target altitude, to constrain end altitude
function constraints = const(dalphadt)
  [StageDynamics] = Stages(dalphadt);
  hend = StageDynamics(end,1);
  global rTarget
%  constraints = [hend-rTarget];
  constraints = [hend-rTarget; StageDynamics(end,4)];
endfunction

% function to define maximum altitude, this helps to guide the rocket into the correct orbital position
function ineq = inequalities(dalphadt)
  [StageDynamics] = Stages(dalphadt);
  global rTarget
  ineq = [-StageDynamics(:,1)+rTarget+1];
  %; StageDynamics(:,5)+15*pi/180; -StageDynamics(:,5)+15*pi/180
endfunction












