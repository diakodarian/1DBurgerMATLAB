% 
% BurgersEqSolver.m 
%         
%
% Author:   Diako Darian
% Date:     10.07.2015
% 
% 
% 
% Purpose    : BurgersEqSolver solves 1D Burgers' equation
%    
%                        u_t + uu_x = nu u_xx
%
% This equation will be solved by using Chebyshev spectral collocation
% method in the case of non-periodic BC, and by a pseudo-spectral 
% Fourier-Galerkin method in the case of a periodic BC.
%
%
% This code solves the above equation for mixed boundary conditions 
% Choose BC = 0 for periodic BC
% Choose BC = 1 for Dirichlet and/or Neumann BC (default)
%
%     --------------oooooo---------------------
%
% Test problem for Fourier-Galerkin method: 
%   
%           Initial condition u(x,0) = u0 = -R sin(x)
%           where R is the Reynolds number
%           Periodic BC 
% 
%          
%
%     --------------oooooo---------------------
%
% Test problem for Chebyshev spectral collocation method: 
%   
%           Initial condition u(x,0) = u0 = 2x 
%           
%           Dirichlet BC u(-1,t) = -2/(1+2t) and u(1,t) = 2/(1+2t)
% 
%           This problem has analytical solution u(x,t) = 2x/(1+2t) 

BC = 1;

%--------------------------------------------------------------------------
%                    Fourier-Galerkin Method
%--------------------------------------------------------------------------
if BC == 0

 % Grid and initial conditions:
  N = 128; nu = 0.01; L = pi; 
  h = 2*pi/N; x = h*(1:N);
  %x = L*(x-pi)/pi;
  %u0 = sin((x-L)/L*pi); 
  %u0 = exp(-4*x.^2);
  u0 = sin(x);
  tmax = 1.0; 
  %PseudoSpectralFourier(N,L,nu,x,u0,tmax); % With integrating factor
  PseudoSpectralFourierGalerkin(N,L,nu,x,u0,tmax)
 
%--------------------------------------------------------------------------
%                  Chebyshev Spectral Collocation Method
%--------------------------------------------------------------------------
elseif BC == 1

  % Grid and initial conditions:
  N = 16; nu = .01; a = -1; b = 1; L = b - a; 
  x = cos(pi*(0:N)/N); Jacobian = 2/L;
  %b = x<=0;
  %u0 = zeros(1,N+1); u0(b) = -sin(pi*x(b)); 
  %u0 = exp(-8*x.^2); 
  u0 = 2*x; 
  alpha = u0(1); beta = u0(end);
  tmax = .5; 
  timeStepping = ['Euler';'RK4  ']; 
  strL1 = 1:5; integrator = 2; % choose a = 1 for Euler or a = 2 for RK4
  BC = ['Dirichlet'; 'Neumann  '];
  strL2 = 1:9; bc = 1; 
  SpectralChebyshevFFT(N,a,b,nu,x,u0,tmax,timeStepping(integrator,strL1), BC(bc,strL2), alpha, beta, Jacobian)
  
end