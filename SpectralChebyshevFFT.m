%  function [tdata,udata] = SpectralChebyshevFFT(N,a,b,nu,x0, u0,tmax, str,
%  BC, gminus, gplus,jacobian)
% 
%  Solve Burgers' eq. u_t + uu_x + nu u_xx = 0 on [-L,L] by
%  
%  Chebyshev collocation method using FFT 
%
%  Dirichlet BC is used: u(-1)= alpha, u(1) = beta.
%
%         
%
% Author:   Diako Darian
% Date:     10.07.2015
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
%
% MMS Test problem for Chebyshev spectral collocation method: 
%   
%           u(x,t) = cos(x)exp(-2 nu t)
%           Initial condition u(x,0) = u0 = cos(x) 
%           
%           Dirichlet BC u(-1,t) = cos(-1)exp(-2 nu t) and 
%           u(1,t) = cos(1)exp(-2 nu t)
%           
%-------------------ooooooooo----------------------------------------------

function [tdata,udata] = SpectralChebyshevFFT(N,a,b,nu,x0, u0,tmax, str, BC, gminus, gplus,jacobian)

% Set up grid and initial data:
    dt = .1/N^2; x = x0; 
    s1 = 'Euler'; s2 = 'RK4  ';
    BC1 = 'Dirichlet'; BC2 = 'Neumann  ';
    xi = .5*b*(1+x) + .5*a*(1-x);
% Time-stepping by leap frog formula:
    clf, drawnow, set(gcf,'renderer','zbuffer')
    nplots = 7; % Number of plots in time
    nplt = floor((tmax/nplots)/dt); nmax = round(tmax/dt);

% Initial condition:
    v = u0; vold = u0; t = 0;
    udata = zeros(nplots+1,N+1); tdata = zeros(1,nplots+1);
    udata(1,:) = v; tdata(1) = t;

%--------------------------------------------------------------------------
%           Time-stepping using Euler method with Dirichlet BC
%--------------------------------------------------------------------------    
    if strcmp(str,s1) == 1 
        for n = 1:nmax
          t = t+dt;
%           uu1 = chebfft(v.^2)';
%           uu2 = chebfft(chebfft(v))';
%           vnew = vold + 2*dt*(-.5*uu1 + nu*uu2);
          du = computeRHS(v,nu,jacobian, BC);
          vnew = vold + 2*dt*du;
          vold = v; v = vnew;
          % Impose Dirichlet BC:
          if strcmp(BC,BC1) == 1
              v(1) = 0; v(end) = 0;
          end
          if mod(n,nplt) == 0 
              i = floor(n/nplt);
              udata(i+1,:) = v; tdata(i+1) = t;
          end
        end
        surf(x,tdata,udata)%waterfall(x,tdata,udata)
        %axis([-1 1 0 tmax -1.1 1.1]),
        xlabel x, ylabel t, zlabel u
    
        
%--------------------------------------------------------------------------
%           Time-stepping using RK4 method with Dirichlet BC
%--------------------------------------------------------------------------
    elseif strcmp(str,s2) == 1
        
        % RK4 parameters:
        a = [1/6 1/3 1/3 1/6]*dt;
        b = [.5 .5 1]*dt;
        
        for n = 1:nmax              % calculation loop
            t = n*dt;
            v1 = v; v0 = v;
            for rk = 1:4
                du = computeRHS(v,nu,jacobian,BC,x,t);
                if rk < 4
                    v = v0; v = v + b(rk)*du;
                end
                v1 = v1 + a(rk)*du;
            end
            v = v1; 
            % Impose Dirichlet BC:
            if strcmp(BC,BC1) == 1
                % General Dirichlet BC:
                v(1) = gplus*exp(-2*nu*t); v(end) = gminus*exp(-2*nu*t);
                % BC for the test problem: 
                %g = 2./(1+2*t);
                %v(1) = g; v(end) = -g;
            end
            if mod(n,nplt) == 0 
                i = floor(n/nplt);
                udata(i+1,:) = v; tdata(i+1) = t;
            end
        end
        clf, drawnow, set(gcf,'renderer','zbuffer');
        [X,T] = meshgrid(x,tdata);
        
        mesh(X,T,udata)
        xlabel x, ylabel t, zlabel u

        % Plot the analytical solution:
        figure(2);
        Z = cos(X).*exp(-2*nu*T);%(2*X)./(1+2*T);
        mesh(X, T, Z)
        xlabel x, ylabel t, zlabel exact   
        
        
        figure(3);
        error = Z-udata; 
        mesh(X,T,error)
        xlabel x, ylabel t, zlabel error 
    end
    
                                

