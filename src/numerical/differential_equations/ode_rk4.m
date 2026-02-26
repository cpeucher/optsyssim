function [t,y] = ode_rk4(odefun,tspan,y0,tstep)
% Basic 4th-order Runge-Kutta scheme
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a basic 4th-order Runge-Kutta scheme for the
% resolution of ordinary differential equations (ODEs). 
% To be used for the sake of simplicity instead of more elaborate
% (adaptive step) matlab solvers such as ode45.
% The call to ode_rk4 is compatible with that to ode45. 
% In particular, in case of resolution of a system of ODEs, the initial
% conditions and the equations are arranged in the form of a column vector.
% 
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% tspan = [0 5];
% y0 = 0;
% [trk,yrk] = ode_rk4(@(t,y) 2*t,tspan,y0,0.1);
% [tm,ym] = ode45(@(t,y) 2*t,tspan,y0);
%
% tstep = 0.01;
% [trk,yrk] = ode_rk4(@(t,y)vdp(t,y,mu),[0 20],[2;0],tstep);
% [tm,ym] = ode45(@(t,y)vdp(t,y,mu),[0 20],[2;0]);
% 
% function dydt = vdp(t,y,mu)
% dydt = [y(2); mu*(1 - y(1)^2)*y(2) - y(1)];
% end
% 
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% odefun            function handle 
%
% tspan             start and end values of the interval over which the
%                       resolution is carried out 
%                       [2-element real vector]
%
% y0                initial values [column vector]
%                       
%                       y0 is a column vector in which y0(i) contains the
%                       initial value of the i-th differential equation.
%
% tstep             fixed step [real scalar]
% 
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% t                 values at which the solutions are estimated
%                       [column vector]              
%
% y                 estimated solution of the ODE [column vector] or of the
%                       system of ODEs [matrix]
%
%                       y is a matrix with length(t) lines and a number of
%                       columns corresponding to the number of equations.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


t = tspan(1):tstep:tspan(2);
% Create time axis

t = t(:);
% Makes it a column vector

nsteps = length(t);

y = zeros(length(y0),nsteps);

y(:,1) = y0;
% Initial conditions 

for istep = 1:nsteps -1

    ti = t(istep);
    yi = y(:,istep);

    k1 = odefun(ti, yi);
    k2 = odefun(ti + tstep/2, yi + k1*tstep/2);
    k3 = odefun(ti + tstep/2, yi + k2*tstep/2);
    k4 = odefun(ti + tstep, yi + k3*tstep);   

    y(:,istep + 1) = yi + tstep/6*(k1 + 2*k2 + 2*k3 + k4);
    
end

y = y.';

end