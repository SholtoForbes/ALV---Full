function dz = rocketDynamics(z,u,phase)

dalphadt = u(1,:);

h = z(1,:);   %Height


v = z(2,:);   %Velocity
m = z(3,:);   %Mass
gamma = z(4,:);
alpha = z(5,:);


xi = z(6,:);
phi = z(7,:);
zeta = z(8,:);
L = 0*ones(1,length(h));

%xi = 0*ones(1,length(h));
%phi = 0*ones(1,length(h));
%zeta = 0*ones(1,length(h));
%L = 0*ones(1,length(h));


if isnan(gamma)
    gamma = 1.5708;
end

FirstStageThrust = dlmread('FirstStageThrust.txt');
SecondStageThrust = dlmread('SecondStageThrust.txt');
ThirdStageThrust = dlmread('ThirdStageThrust.txt');

switch phase
case 'prepitch'
  T = interp1(FirstStageThrust(:,1),FirstStageThrust(:,6),h/1000)*1000*3*4;
  case 'postpitch'
  T = interp1(FirstStageThrust(:,1),FirstStageThrust(:,6),h/1000)*1000*3*4;
  case 'secondstage'
  T = interp1(SecondStageThrust(:,1),SecondStageThrust(:,6),h/1000)*1000*4;
  case 'thirdstage'
  T = interp1(ThirdStageThrust(:,1),ThirdStageThrust(:,6),h/1000)*1000*2;
end


density = 1.474085291*(0.9998541833.^h);  %Data fit off of wolfram alpha

speedOfSound = 280;  %(m/s)  %At 10 km altitude
mach = v/speedOfSound;

switch phase
case 'prepitch'
  A = .28 + 4*0.5; % Reference Area of first stage with 4 boosters, each booster is 0.5 and core stage is 0.28 (m^2);
  case 'postpitch'
  A = .28 + 4*0.5; % Reference Area of first stage with 4 boosters, each booster is 0.5 and core stage is 0.28 (m^2);
  case 'secondstage'
  A = 0.28;
  case 'thirdstage'
  A = 0.28;
end

%%%% Compute the drag:
AeroCoeffs = dlmread('AeroCoeffs.txt');

Cd = interp1(AeroCoeffs(:,1),AeroCoeffs(:,2),mach);
D = 0.5*Cd.*A.*density.*v.^2;

%%%% Compute gravity from inverse-square law:
rEarth = 6.3674447e6;  %(m) radius of earth
mEarth = 5.9721986e24;  %(kg) mass of earth
G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
g = G*mEarth./((h+rEarth).^2);

switch phase
  case 'prepitch'
  dm = -16.39*4;
  case 'postpitch'
  dm = -16.39*4;
  case 'secondstage'
  dm = -3.952;
  case 'thirdstage'
  dm = -0.4744;
end



switch phase
    case 'prepitch'
    gamma = 1.5708*ones(1,length(h)); % Control Trajectory Angle 
    case 'postpitch'
    %Do nothing
    case 'secondstage'
    %Do nothing
    case 'thirdstage'
    %Do nothing
end


[dr,dxi,dphi,dgamma,dv,dzeta] = RotCoords(h+rEarth,xi,phi,gamma,v,zeta,L,D,T,m,alpha,phase);

switch phase
    case 'prepitch'
    dz = [dr;dv;dm;0;dalphadt;dxi;dphi;0];
    case 'postpitch'
    dz = [dr;dv;dm;dgamma;dalphadt;dxi;dphi;dzeta];
    case 'secondstage'
    dz = [dr;dv;dm;dgamma;dalphadt;dxi;dphi;dzeta];
    if h < 0 % This section limits the simulation from straying into computationally bad territory
    dz = zeros(8,1);
    elseif h > 1000000
    dz = zeros(8,1);
    elseif gamma < -1.57
    dz = zeros(8,1);
    elseif gamma > 1.57
    dz = zeros(8,1);
    end
    case 'thirdstage'
    dz = [dr;dv;dm;dgamma;dalphadt;dxi;dphi;dzeta];
    if h < 0
      dz = zeros(8,1);
    elseif h > 1000000
      dz = zeros(8,1);
    elseif gamma < -1.57
    dz = zeros(8,1);
    elseif gamma > 1.57
    dz = zeros(8,1);
    end
    
    
    
    
end


end