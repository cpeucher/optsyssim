% -------------------------------------------------------------------------
% Test of ode_rk4 function
%
% We conduct simple standard tests for the resolution of ordinary
% differential equations (ODEs) and systems of ODEs using the basic 4th
% order Runge-Kutta (RK4) scheme, and compare the results with those 
% obtained with the Matlab ode45 solver (default parameters).
%
% Play with the step-size for the RK4 scheme.
%
% Test 1: quadratic function
% Test 2: van der Pol oscillator
% Test 3: decaying exponential
% Test 4: Lorenz attractor
% Test 5: sinusoidal function
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all
format long

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
time_stamp = string(datetime('now','TimeZone','Z','Format','yyyyMMdd''T''HHmmss''Z'));

% ------------------------------------------------------------------------- 
% Reinitialise the random number generator for reproducibility.
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
reset(stream);

%--------------------------------------------------------------------------
% Switches
%--------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
do_add_figsize_to_filename = 1;
margin_figure = 0;

% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
fig.interpreter = 'latex';
fig.font_size = 18;

mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));



% -------------------------------------------------------------------------
% Test 1: quadratic function
% -------------------------------------------------------------------------

tspan = [0 5];
y0 = 0;

fprintf('\n%s','Test 1: quadratic function')
fprintf('\n%s','==========================')

fprintf('\n\n%s\n','Matlab ode45')

tic
[tm,ym] = ode45(@(t,y) 2*t,tspan,y0);
toc
% Matlab solver

em = rmse(ym,tm.^2);

fprintf('%s%3.2e','RMSE: ',em)


fprintf('\n\n%s\n','ode_rk4')

tic
[trk,yrk] = ode_rk4(@(t,y) 2*t,tspan,y0,0.1);
toc
% 4-th order Runge-Kutta scheme

figure('Name','Test 1')
plot(tm,tm.^2,'b-')
hold on
plot(tm,ym,'ro')
plot(trk,yrk,'k^')
legend('exact','ode45','ode\_rk4')
xlabel('t')
ylabel('y(t)')

errk = rmse(yrk,trk.^2);

fprintf('%s%3.2e\n','RMSE: ',errk)


%%
% -------------------------------------------------------------------------
% Test 2: van der Pol ODE
% -------------------------------------------------------------------------

mu = 1;
% Parameter of van der Pol ODE

tstep = 0.01;
% Step-size for RK4 scheme


fprintf('\n%s','Test 2: van der Pol ODE')
fprintf('\n%s','=======================')

fprintf('\n\n%s\n','Matlab ode45')
tic
[tm,ym] = ode45(@(t,y)vdp(t,y,mu),[0 20],[2;0]);
toc

fprintf('\n\n%s\n','ode_rk4')

tic
[trk,yrk] = ode_rk4(@(t,y)vdp(t,y,mu),[0 20],[2;0],tstep);
toc
% 4-th order Runge-Kutta scheme

figure('Name','Test 2')
plot(tm,ym(:,1),'-ob')
hold on
plot(tm,ym(:,2),'-or')
plot(trk,yrk(:,1),'*b')
plot(trk,yrk(:,2),'*r')
legend('ode45 - y1','ode45 - y2','ode\_rk4 - y1','ode\_rk4 - y2')
xlabel('t')
ylabel('y(t)')


function dydt = vdp(t,y,mu)
% van der Pol ODE

dydt = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];

end

%%
% -------------------------------------------------------------------------
% Test 3: somehow stiffer problem - decaying exponential
% -------------------------------------------------------------------------

tspan = [0 0.02];
y0 = 1;

fprintf('\n%s','Test 3: decaying exponential')
fprintf('\n%s','============================')

fprintf('\n\n%s\n','Matlab ode45')

tic
[tm,ym] = ode45(@test3,tspan,y0);
toc
% Matlab solver

em = rmse(ym,exp(-1000*tm));

fprintf('%s%3.2e','RMSE: ',em)


fprintf('\n\n%s\n','ode_rk4')

tic
[trk,yrk] = ode_rk4(@test3,tspan,y0,0.001);
toc
% 4-th order Runge-Kutta scheme

figure('Name','Test 3')
plot(tm,exp(-1000*tm),'b-')
hold on
plot(tm,ym,'ro')
plot(trk,yrk,'k^')
legend('exact','ode45','ode\_rk4')
xlabel('t')
ylabel('y(t)')

figure('Name','Test 3 - log')
semilogy(tm,exp(-1000*tm),'b-')
hold on
semilogy(tm,ym,'ro')
semilogy(trk,yrk,'k^')
legend('exact','ode45','ode\_rk4')
xlabel('t')
ylabel('y(t)')

errk = rmse(yrk,exp(-1000*trk));

fprintf('%s%3.2e\n','RMSE: ',errk)


function dydt = test3(t,y)
dydt = -1000*y;
end


%%
% -------------------------------------------------------------------------
% Test 4: Lorenz attractor
% -------------------------------------------------------------------------

tspan = [0 100];
y0 = [1 1 1];

fprintf('\n%s','Test 4: Lorenz attractor')
fprintf('\n%s','========================')

fprintf('\n\n%s\n','Matlab ode45')

tic
[tm,ym] = ode45(@lorenz,tspan,y0);
toc
% Matlab solver



fprintf('\n\n%s\n','ode_rk4')

tic
[trk,yrk] = ode_rk4(@lorenz,tspan,y0,0.001);
toc
% 4-th order Runge-Kutta scheme


figure('Name','Test 4')
plot3(ym(:,1),ym(:,2),ym(:,3),'b-')
hold on
plot3(yrk(:,1),yrk(:,2),yrk(:,3),'r--')
legend('ode45','ode\_rk4')
xlabel('x')
ylabel('y')
zlabel('z')




function dydt = lorenz(t,y)
% Lorenz system

sigma = 10;
rho = 28;
beta = 8/3;

dydt(1) = sigma*(y(2) - y(1));
dydt(2) = y(1)*(rho - y(3)) - y(2);
dydt(3) = y(1)*y(2) - beta*y(3);

dydt = dydt(:);

end


%%
% -------------------------------------------------------------------------
% Test 5: sinusoidal function
% -------------------------------------------------------------------------

tspan = [0 3*pi];
y0 = 0;

fprintf('\n%s','Test 5: sinusoidal function')
fprintf('\n%s','===========================')

fprintf('\n\n%s\n','Matlab ode45')

tic
[tm,ym] = ode45(@(t,y) cos(t),tspan,y0);
toc
% Matlab solver

em = rmse(ym,sin(tm));

fprintf('%s%3.2e','RMSE: ',em)


fprintf('\n\n%s\n','ode_rk4')

tic
[trk,yrk] = ode_rk4(@(t,y) cos(t),tspan,y0,0.1);
toc
% 4-th order Runge-Kutta scheme



errk = rmse(yrk,sin(trk));

fprintf('%s%3.2e\n','RMSE: ',errk)

figure('Name','Test 5')
plot(tm,sin(tm),'b-')
hold on
plot(tm,ym,'ro')
plot(trk,yrk,'k^')
legend('exact','ode45','ode\_rk4')
xlabel('t')
ylabel('y(t)')










% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
% alignfigs(2)


% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

