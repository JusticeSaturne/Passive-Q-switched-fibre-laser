function [Psf,time,N2,Nsa2] = PassiveQswitch(Pump,Leng,Nt)
% Function PassiveQswitch implementes a finite difference method to compute
% output characteristics of an Ytterbiium doped passive Q-switched fibre laser
% using Cr4+:YAG saturable absorber
%
% Inputs:
%  - Pump: Input pump power
%  - Leng:Length of the Ytterbium doped fibre
%  - Nt: doped fibre concentration
%  
%
% Outputs:
%   
%   - N2: excited state population density variation with time
%   - Nsa2 : excited state 
%   - Psf: Q-switched fibre laser output power
%   - time: computation time
% Comments:
%   - Ytterbium doped gain medium is a two level system
%
% References:
%
% Written by Justice Sompo, University of Johannesburg, South Africa
%
%
%
close all
clc
% Check inputs
if nargin < 3
    Nt = 9e25;
    if nargin < 2
        Leng = 5;
        if nargin < 1
            error(message('MATLAB:PassiveQswitch:NotEnoughInputs'));
        end
    end
end
lambdas = 1064e-9;              % Laser signal waveength
sigmaas = 1.16e-27;             % Absorption cross section at signal wavelength
sigmaes = 2.34e-25;             % Emission cross section at signal wavelength
sigmagsa = 3.74e-22;            % Absorption cross section of the saturable absorber 
sigmaesa = 0.935e-22;           % Emission cross section of the saturable absorber
NA = 0.15;
dcore = 5.4e-6;                 % fiber core diameter
dclad = 125e-6;                 % Fiber cladding diameter                
Step_number = 50;               % Number of section
L = Leng;                       % Length of the doped fibre 
R = 0.04;
alfap = 0.005;                  % Background loss at pump wavelength
alfas = 0.005;                  % Background loss at laser wavelength
alfasa = 20;                    % Background loss of the saturable absorber
nfiber = 1.45;                  % Refractive index of the fibre
n0 = Nt;                        % Concentration of doped fibre
c0 = 3e8;
V = 2*pi*dcore*NA/(2*lambdas);
w = (dcore/2)*(0.65+1.619*V^(-1.5)+2.879*V^(-6));
gamas = 1-exp(-(2*dcore^2/(4*w^2)));
wsa0 = 7.2e-6;                   % Focused beam waist in SA
nsa = 6.22e23;
nSA = 1.8;
taosa = 1e-6;
lsa = 0.0025;
Step_number1 = 25;
z0 = pi*wsa0*wsa0*nSA/(lambdas);
%==========================================================================
q = 1:Step_number1+1;
z(q) = (lsa/(Step_number1-1))*(q-(2+Step_number1)/2);
wsa(q) = wsa0*sqrt(1+(z(q)/z0).^2);
Asa(q) = pi*wsa(q).^2;
%==========================================================================
h = 6.626e-34;                 %Plank constant
vg = c0/nfiber;
deltax1 = lsa/Step_number1;
deltax = L/Step_number;
deltat = L/(Step_number*vg);
deltalamda = 3e-6;
Aco = pi*(dcore/2)^2;        % Doped fibre core diameter
lamdap = 976e-9;             % Pump wavelength
gamap = dcore^2/(dclad^2);   % Overlap factor between the pump and doped area
sigmaep = 2.44e-24;          % Emission cross section at pump power
sigmaap = 2.5e-24;           % Absorption cross section at pump power
tao = 850e-6;
%==========================================================================
% Initial conditions of the saturable absorber
ssa = 1:Step_number1+1;
nsa2(ssa,1) = 100;
psfSA(ssa,1) = 0;
psbSA(ssa,1) = 0;
%==========================================================================
% Initial condition for the Ytterbum doped fibre 
s = 1:Step_number+1;
n2(s,1) = 100;
ppb(s,1) = 0.0;
psf(s,1) = 0.0;
psb(s,1) = 0.0;
%==========================================================================

Pp = Pump;              %(Pump power in watt)
temps = deltat*450000;  % Simulation time value can be choosen arbitrary
Time = deltat:deltat:temps;
ppb(Step_number,1) = Pp;
t = 0;
for ii = 1:length(Time)
    ppb(Step_number+1,2) = Pp;
    for s = 1:Step_number
        n2(s,2) = n2(s,1)+deltat*((gamap*lamdap/(h*c0*Aco))*(sigmaap*(n0-n2(s,1))...
            -sigmaep*n2(s,1))*ppb(s,1)+(gamas*lambdas/(h*c0*Aco))*(sigmaas*(n0-n2(s,1))...
            -sigmaes*n2(s,1))*(psf(s,1)+psb(s,1))-n2(s,1)/tao);
        ppb(s,2) = ppb(s+1,1)+deltax*(gamap*(sigmaep*n2(s+1,1)...
            -sigmaap*(n0-n2(s+1,1)))*ppb(s+1,1)-alfap*ppb(s+1,1));
        psb(s,2) = psb(s+1,1)+deltax*(gamas*(sigmaes*n2(s+1,1)...
            -sigmaas*(n0-n2(s+1,1)))*psb(s+1,1)+2*gamas*sigmaes*n2(s+1,1)...
            *h*c0*c0*deltalamda/(lambdas)^3-alfas*psb(s+1,1));
        psf(s+1,2) = psf(s,1)+deltax*(gamas*(sigmaes*n2(s,1)-sigmaas*(n0-n2(s,1)))...
            *psf(s,1)+2*gamas*sigmaes*n2(s,1)*h*c0*c0*deltalamda/(lambdas)^3-alfas*psf(s,1));
    end
    n2(Step_number+1,2) = n2(Step_number+1,1)+deltat*((gamap*lamdap/(h*c0*Aco))...
        *(sigmaap*(n0-n2(Step_number+1,1))-sigmaep*n2(Step_number+1,1))*ppb(Step_number+1,1)...
        +(gamas*lambdas/(h*c0*Aco))*(sigmaas*(n0-n2(Step_number+1,1))...
        -sigmaes*n2(Step_number+1,1))*(psf(Step_number+1,1)+psb(Step_number+1,1))...
        -n2(Step_number+1,1)/tao);
    for s1=1:Step_number1
        nsa2(s1,2) = nsa2(s1,1)+deltat*(((psfSA(s1,1)+psbSA(s1,1))*sigmagsa*lambdas)...
            *(nsa-nsa2(s1,1))/(Asa(s1)*h*c0)-nsa2(s1,1)/taosa);
        psbSA(s1,2) = psbSA(s1+1,1)-deltax1*(sigmaesa*nsa2(s1+1,1)...
            +sigmagsa*(nsa-nsa2(s1+1,1))-alfasa)*psbSA(s1+1,1);
        psfSA(s1+1,2) = psfSA(s1,1)-deltax1*(sigmaesa*nsa2(s1,1)...
            +sigmagsa*(nsa-nsa2(s1,1))-alfasa)*psfSA(s1,1);
    end
    nsa2(Step_number1+1,2) = nsa2(Step_number1+1,1)+deltat*(((psfSA(Step_number1+1,1)...
        +psbSA(Step_number1+1,1))*sigmagsa*lambdas)*(nsa-nsa2(Step_number1+1,1))/(Asa(Step_number1+1)...
        *h*c0)-nsa2(Step_number1+1,1)/taosa);
    psbSA(Step_number1+1,2) = 0.85*psb(1,1);
    psfSA(1,2) = psbSA(1,1);
    psf(1,2) = 0.85*psfSA(Step_number1+1,1);
    psb(Step_number+1,2) = psf(Step_number+1,1)*R*0.85;
    for s1 = 1:Step_number1+1
        nsa2(s1,1) = nsa2(s1,2);
        psbSA(s1,1) = psbSA(s1,2);
        psfSA(s1,1) = psfSA(s1,2);
    end
    for s = 1:Step_number+1
        n2(s,1) = n2(s,2);
        ppb(s,1) = ppb(s,2);
        psb(s,1) = psb(s,2);
        psf(s,1) = psf(s,2);
    end
    t = t + deltat;
    Psf(ii) = psf(Step_number+1);
    Psb(ii) = psb(Step_number+1);
    Ppb(ii) = ppb(Step_number+1);
    time(ii) = t;
    N2(ii) = n2(Step_number+1);
    Nsa2(ii) = nsa2(Step_number1+1);
end


