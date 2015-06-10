% function [tdata,udata] = PseudoSpectralFourier(N,L,nu,x0, u0,tmax)
% Solve Burgers' eq. u_t + uu_x + nu u_xx = 0 on [-L,L] by
% FFT with integrating factor v = exp(-nu k^2t)*u-hat.
% 
% PseudoSpectralFourier.m 
%         
%
% Author:   Diako Darian
% Date:     10.07.2015
% 
% 
% 
% Purpose    : PseudoSpectralFourier.m solves 1D Burgers' equation
%    
%                        u_t + uu_x = nu u_xx
%
% by FFT with integrating factor v = exp(-nu k^2t)*u-hat.
%
%
%     --------------oooooo---------------------
%
% Test problem for Fourier-Galerkin method: 
%   
%   
%           Initial condition u(x,0) = u0 = -R sin(x),
%           where R is the Reynolds number
%           Periodic BC 
%
%     --------------oooooo---------------------

function [tdata,udata] = PseudoSpectralFourier(N,L,nu,x0, u0,tmax)

% Set up grid and initial data:
  dt = 0.0005;%.1/N^2;
  h = 2*pi/N; x = x0;
  
% Initial condition:
  clf, drawnow, set(gcf,'renderer','zbuffer')
  u = u0;
  v = fft(u); k = [0:N/2-1 0 -N/2+1:-1]*(pi/L); k2 = k.^2;

% Filter for dealiasing nonlineat convection:
  kmax = (2/3)*(N/2 + 1)*(pi/L);
  dealias = abs(k)>kmax;

% Solve PDE and plot results:
  nplt = floor((tmax/7)/dt); nmax = round(tmax/dt);
  udata = zeros(8,N); tdata = zeros(1, 8);
  udata(1,:) = u; tdata(1) = 0; h = waitbar(0,'please wait...');
  for n = 1:nmax
    t = n*dt; g = -.5i*dt*k;
    E = exp(-dt*nu*k2/2); E2 = E.^2;
    a = g.*fft(real( ifft(     v    ) ).^2);
    b = g.*fft(real( ifft(E.*(v+a/2)) ).^2);     % 4th-order
    c = g.*fft(real( ifft(E.*v + b/2) ).^2);     % Runge-Kutta
    d = g.*fft(real( ifft(E2.*v+E.*c) ).^2);
    v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
    v(dealias) = 0;
    if mod(n,nplt) == 0 
      u = real(ifft(v)); waitbar(n/nmax)
      i = floor(n/nplt);
      udata(i+1,:) = u; tdata(i+1) = t;
    end
  end
  waterfall(x,tdata,udata), xlabel x, ylabel t, zlabel u
  %waterfall(x,tdata,udata), colormap default; view(-20,25)
  xlabel x, ylabel t, zlabel u%, axis([-L+h L 0 tmax -1.1 1.1]), grid off
%   udata(1, N/2)
%   udata(end, N/2)
%   figure;
%   [X,T] = meshgrid(x,tdata);
%   Z = (2*X)./(1+2*T);
%   surf(X, T, Z)
%   xlabel x, ylabel t, zlabel Z
%   error = zeros(8,N); 
%   for i = 1:8
%     error(i,:) = abs(Z(i,:)-udata(i,:));
%   end
%   figure(3);
%   surf(x,tdata,error)
