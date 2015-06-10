%  function [getRHS] = computeRHS(v, nu, jacobian, BC)
% 
%
%  This function computes the RHS of Burger's equation
% 
%  Compute RHS: du = -0.5(u^2)_x + nu u_xx 
%
%         
%
% Author:   Diako Darian
% Date:     10.07.2015
% 
% 
% 
% Purpose    : computeRHS.m comutes RHS of 1D Burgers' equation
%    
%                        u_t + uu_x = nu u_xx
%
% by FFT.
%
%-------------------ooooooooo----------------------------------------------

function [getRHS] = computeRHS(v, nu, jacobian, BC)

BC1 = 'Dirichlet'; BC2 = 'Neumann  ';
gplus = 0; gminus = 0;
% Diffusion term:
% Impose Neumann BC:
if strcmp(BC,BC2) == 1
    uu2 = chebfft(v);
    uu2(1) = gplus; uu2(end) = gminus;
    u2 = chebfft(uu2)';
elseif strcmp(BC,BC1) == 1
    u2 = chebfft(chebfft(v))';
end


%--------------------------------------------------------------------------
% Nonlinear term:
% Two approaches:

% 1)  This eqution gives aliasing error:
u1 = chebfft(v.^2)';

%..........................................................................

% 2)  The alternative approach for the non-linear term
%     that removes the aliasing error
N = length(v)-1;
x = cos((0:N)'*pi/N);
ii = 0:N-1;
vv = v.^2; vv = vv(:); 
v_hat = fft(vv); 

% Filter for dealiasing nonlineat convection:
k = [0:N/2-1 0 -N/2+1:-1];
kmax = (2/3)*(N/2 + 1);
dealias = abs(k)>kmax;
v_hat(dealias) = 0;
vr = real(ifft(v_hat));

V = [vr; flipud(vr(2:N))];      % transform x -> theta          
U = real(fft(V));
W = real(ifft(1i*[ii 0 1-N:-1]'.*U));
w = zeros(N+1,1);
w(2:N) = -W(2:N)./sqrt(1-x(2:N).^2);    % transform theta -> x     
w(1) = sum(ii'.^2.*U(ii+1))/N + .5*N*U(N+1);     
w(N+1) = sum((-1).^(ii+1)'.*ii'.^2.*U(ii+1))/N + ...
          .5*(-1)^(N+1)*N*U(N+1);
%--------------------------------------------------------------------------

% assemble the RHS:
getRHS = -.5*jacobian*w' + nu*(jacobian^2)*u2;  
%getRHS = -.5*u1 + nu*u2;