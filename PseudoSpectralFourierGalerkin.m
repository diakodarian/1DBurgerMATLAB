%
%function [tdata,udata] = PseudoSpectralFourierGalerkin(N,L,nu,x0,v0,tmax)
%
% Solve Burgers' eq. u_t + uu_x + nu u_xx = 0 on [-L,L] by FFT.
% 
%
% PseudoSpectralFourierGalerkin.m 
%         
%
% Author:   Diako Darian
% Date:     14.07.2015
% 
% 
% 
% Purpose    : PseudoSpectralFourierGalerkin.m solves 1D Burgers' equation
%    
%                        u_t + uu_x = nu u_xx
%
% by FFT.
%
%
%     --------------oooooo---------------------
%
% MMS Test problem for Fourier-Galerkin method: 
%   
%   
%           Initial condition u(x,0) = u0 = sin(x),
%           Analytic solution u(x,t) = sin(x)exp(-2 nu t),
%           Periodic BC 
%
%     --------------oooooo---------------------



%------------------------ooooooooooooo-------------------------------------
%                  The Fourier-Galerkin Approach
%--------------------------------------------------------------------------


function [tdata,udata] = PseudoSpectralFourierGalerkin(N,L,nu,x0,v0,tmax)

% Set up grid and initial data:
  dt =0.0001;%0.5/N^2;
  x = x0;
  
% Initial condition:
  u = fft(v0); 

% Wave number vector:  
  k = [0:N/2-1 0 -N/2+1:-1]*(pi/L); k2 = k.^2;
  
% Filter for dealiasing nonlineat convection:
  kmax = (2/3)*(N/2 + 1)*(pi/L);
  dealias = abs(k)>kmax;

%--------------------------------------------------------------------------
%           Time-stepping using RK4 method
%--------------------------------------------------------------------------

    % RK4 parameters:
    a_rk4 = [1/6 1/3 1/3 1/6]*dt;
    b_rk4 = [.5 .5 1]*dt; 
    
    nplot = 20; % number of plots in time
    clf, plotgap = round(tmax/dt); nplt = floor((tmax/nplot)/dt);
    tdata = zeros([1 nplot]); udata = zeros([nplot, N]);
    tdata(1, 1) = 0; udata(1,:) = v0;
    u_exact = zeros([nplot, N]); u_exact(1,:) = v0; i = 0;
    for n = 1:plotgap
        t = n*dt; 
        % RK4 loop 
        u1 = u; u0 = u;
        for rk = 1:4
            du = RHS_Burger(nu,k,k2,dealias,x,u,t);
            if rk < 4
                u = u0; u = u + b_rk4(rk)*du;
            end
            u1 = u1 + a_rk4(rk)*du;
        end % end for RK4 loop
        u = u1;
        % Store data
        if mod(n,nplt) == 0 
            tdata(1, i+1) = t; udata(i+1,:) = real(ifft(u)); 
            u_exact(i+1, :) = sin(x).*exp(-2*nu*t);
            i = i+1;
        end
        
    end % end for (time integration)

 %             -------------oooooooo----------------  

%------------------------ooooooooooooo-------------------------------------
%                           MMS Test
%--------------------------------------------------------------------------

clf, drawnow, set(gcf,'renderer','zbuffer');

[X,T] = meshgrid(x,tdata);
% mesh(X,T,u_exact);
% xlabel x, ylabel t, zlabel u_exact
%figure(2);
%mesh(X,T,udata);
%xlabel x, ylabel t, zlabel u

%figure(3);
err = udata-u_exact;
mesh(X,T, err);
xlabel x, ylabel t, zlabel error
%------------------------ooooooooooooo-------------------------------------


%--------------------------------------------------------------------------
%           RHS of 1D Burgers' equation
%--------------------------------------------------------------------------

function [rhs] = RHS_Burger(nu,k,k2,dealias,x,u_hat,t)

    %----------------------------------------------------------------------
    %                     Non-linear term
    u_real = real(ifft(u_hat)); u_real2 = u_real.^2; u_real2_hat = fft(u_real2);
    nonlinear_term = 1i*k.*u_real2_hat;
    nonlinear_term(dealias) = 0;
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %                     Viscous term
    viscous_term = nu*k2.*u_hat;
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    %                     MMS term
    f = -nu*sin(x)*exp(-2*nu*t) + sin(x).*cos(x)*exp(-4*nu*t);
    f_hat = fft(f);
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    %                     Assamble RHS
    rhs = -.5*nonlinear_term - viscous_term + f_hat;
    %----------------------------------------------------------------------
    
    





