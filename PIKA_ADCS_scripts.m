%% PIKA Disturbance Torques Script

%The functions below derived from: 
%J. Eterno and S. Starin, NASA, tech.,
%Jan. 2011.https://ntrs.nasa.gov/api/citations/20110007876/downloads/20110007876.pdf.

%Note: Most of the actuator-sizing relevant functions below are for
%reaction wheels. Let me know if support for other actuator functions are needed. For
%more information on both disturbance torques and actuator sizing see pgs.
%8-11 and 18-23 of https://ntrs.nasa.gov/api/citations/20110007876/downloads/20110007876.pdf

%List of Functions:
% --------------------------------------------------------
    %SRP --> estimates the solar radiation pressure disturbance torque given 6
    %parameters. The parameters are: r,h,psi,dc,q,phi
        %r [m] --> radius of the spacecraft (assuming it is ~cylindrical)
        %h [m]--> height of the spaceraft
        %psi [deg] -->angle of incidence, if you let psi == 0, the function will use a
        %default value
        %dc [m] --> distance between center of solar radiation pressure and center of
        %mass, if you let dc == 0, the function will use a default value
        %q [unitless] --> reflectance factor --> if you let q ==0, the function will use a
        %default value
        %phi [W/m^2] --> solar constant --> also has a default value for an input of phi
        %==0
        %NOTE: the defaulted values assume a small spacecraft in an LLO orbit, so
        %if you want to calculate a disturbance torque for s/c in a different
        %orbital trajectory or of a different size, make sure to adjust the
        %parameter values.
        
    %Grav --> estimates the gravity gradient torque given 6 paramaters. The
    %parameters are mu,R,theta,r,h,m
        %mu [m^3s^-2] --> gravitational parameter. If you let mu == 0,
        %it will default to mu for the moon
        %R [m] --> distance from orbit to center of orbiting body. If you let
        %R == 0 , it will default to R given an orbit with altitude 100*10^3 m
        %and the orbiting body is the moon
        %theta [deg] --> distance between the local vertical and the Z
        %principal axis. if you let **theta == -1**, it will default to 10
        %degrees.
        %r,h are as defined above
        %m [kg] --> is the mass of the spacecraft
        %NOTE: the defaulted values assume a small spacecraft in an LLO orbit, so
        %if you want to calculate a disturbance torque for s/c in a different
        %orbital trajectory or of a different size, make sure to adjust the
        %parameter values.
        
    %Control --> estimates a control torque given a disturbance torque and
    %a safety margin. The parameters are Td, margin
        %Td [N*m]--> the disturbance torque from which to calculate a control 
        %margin [unitless] --> safety margin. If margin == 0, will default to safety
        %margin of 1.2.
        
    
    %Slew --> estimates the torque required for a certain slew. Assumes max
    %acceleration slews with no resisting momentum (1/2 distance in 1/2
    %time). The parameters are theta, I, t
        %theta [deg] --> desired slew (i.e. 30 deg slew, 10 deg slew,
        %etc...)
        %I [kg*m^2] --> moment of inertia for the axis about which to make
        %the desired slew
        %t [s] --> desired slew time
    
     %Moment --> estimates the maximum momentum saturation given a certain
     %orbit. Assumes max saturation 1/4 of the way through the orbit.
     %Parameters are Td, P
        %Td [N*m] --> disturbance torque
        %P [s] --> period of the orbit

    
%% Testing functions/Examples

TS = SRP(1.4,0.9,0,0,0,0); %spacecraft with r =1.4m and h = 0.9 m and defaulted values otherwise

TG = Grav(0,0,-1,1.4,0.9,70); %same dimensions as above and mass = 70 kg

Tc = Control(TG,0); %uses the primary disturbance torque in this example, which is gravity gradient torque

h = Moment(TG,0); %finds an estimate of the maximum reaction wheel saturation



%% functions start here

%solar radiation pressure disturbance torque
function Ts = SRP(r,h,psi,dc,q,phi)%please input SI units
c = 3*10^8; %speed of light [m/s]

%The value of phi, the Total Solar Irradiance (TSI) or solar constant on
%the moon varies depending on the location of moon relative to the earth.
%The average TSI below just finds the mean of all the possible extrema.
%See: %https://www.e3s-conferences.org/articles/e3sconf/pdf/2018/24/e3sconf_solina2018_00053.pdf
%for more details
TSI = mean([1400.83 1416.36 1311.18 1323.74 1310.44 1324.49 1401.65 1415.54]);

%default to mean TSI
if phi == 0
    phi = TSI;
end

%default to angle of incidence of 0
if psi == 0
    psi = 0;
end

%default to reflectance of 0.6
if q == 0
    q = 0.6;
end

%default to distance between center of radiation pressure and center of
%mass of 0.1. I chose a small dc as the default because our s/c is
%relatively small. 
if dc == 0
    dc = 0.1
end
%Assume the maximum sunlit surface area is a semicylinder
As = pi*r*h + pi*r^2 + 2*r*h;


Ts = (phi/c)*As*(1+q)*(dc)*cos(psi);

end


%gravity gradient disturbance torque
function Tg = Grav(mu,R,theta,r,h,m)

rm = 1738*10^3;

%Gravitational parameter --> defaults to mu for the moon
if mu == 0
    mu = 4.905*10^12;
end

%distance from the center of the orbiting body --> default assumes the
%orbiting body is the moon with an orbital altitude of 100 km
if R == 0
    R = 100*10^3 + rm;
end

%angle between local vertical and Z principal axis --> assuming a larger
%difference between geometric and principal axes to accomodate weight
%distrubution imbalances due to tertary payloads, ESPA ring, etc..
if theta == -1
    theta = 10;
end

Iz = 0.5*m*r^2; 
Iy = (1/12)*m*(3*r^2 + h^2);

Tg = ((3*mu)/(2*R^3))*abs(Iz - Iy)*sind(2*theta);

end


%Control torque
function Tc = Control(Td,margin)

%default safety margin of 1.2 --> just an estimate
if margin == 0
    margin = 1.2
end


Tc = Td*margin;
end

%Slew momentum
function slw = Slew(theta,I,t)
slw = (4*theta*I)/(t^2);
end


%Max Momentum saturation
function h = Moment(Td,P)

if P == 0
    P = 2*3600
end

%rms average of sinusoidal function divided by 4 (assume max saturation 1/4 way
%through orbit)
c = 0.707/4;

h = Td*P*c;
end





