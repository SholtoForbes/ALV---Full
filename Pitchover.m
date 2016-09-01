%Pitchover angle Computation
function pitchover_angle = Pitchover(icond,vehicle,rTarget_FirstStage,prepitch_time)

h0_prepitch = icond.h0_prepitch;
xi0_prepitch = icond.xi0_prepitch;
phi0_prepitch = icond.phi0_prepitch;
gamma0_prepitch = icond.gamma0_prepitch;
zeta0_prepitch = icond.zeta0_prepitch;

   
mFirstStage = vehicle.mFirstStage;
mFirstStageFuel = vehicle.mFirstStageFuel;

mSecondStage = vehicle.mSecondStage;
mSecondStageFuel = vehicle.mSecondStageFuel;
mThirdStage = vehicle.mThirdStage;
mThirdStageFuel = vehicle.mThirdStageFuel;
mPayload = vehicle.mPayload;

mTotal = vehicle.mTotal;

mdotFirstStage = vehicle.mdotFirstStage;
mdotSecondStage = vehicle.mdotSecondStage;
mdotThirdStage = vehicle.mdotThirdStage;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Pre-Pitchover Simulation                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
printf('Calculating Pitchover Angle    ')
fflush(stdout);

printf('Pre-Pitch Running    ')
fflush(stdout);

v0_prepitch = 0;  %Rocket starts stationary
m0_prepitch = mTotal;  %Rocket starts full of fuel

phase = 'prepitch';
tspan1 = [0:prepitch_time]; % prepitch time 
prepitch0 = [h0_prepitch, v0_prepitch, m0_prepitch, gamma0_prepitch, 0, xi0_prepitch, phi0_prepitch, zeta0_prepitch];

[prepitch] = lsode(@(prepitch,t) rocketDynamics(prepitch,t,0,phase,tspan1), prepitch0, tspan1);  


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Post-Pitchover Simulation                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
printf('Post-Pitch Running    ')
fflush(stdout);

pitchover_angle = fminbnd(@(x) postpitch(x,rTarget_FirstStage,mFirstStageFuel,mdotFirstStage,prepitch),pi/180*0.1,pi/180*10)

printf('Complete    ')
fflush(stdout);

endfunction

function err = postpitch(pitchover_angle,rTarget_FirstStage,mFirstStageFuel,mdotFirstStage,prepitch)

  phase = 'postpitch';
  tspan2 = [15:mFirstStageFuel*4/(mdotFirstStage*4)]; % run until fuel ends
  tspan2(end) = mFirstStageFuel*4/(mdotFirstStage*4);

  postpitch0 = [prepitch(end,1), prepitch(end,2), prepitch(end,3), pi/2-pitchover_angle, 0, prepitch(end,6), prepitch(end,7), prepitch(end,8)];
  dalphadt = 0; % gravity turn, AoA=0

  [postpitch] = lsode(@(postpitch,t) rocketDynamics(postpitch,t,dalphadt,phase,tspan2), postpitch0, tspan2);

  err = abs(postpitch(end,1) - rTarget_FirstStage);
endfunction