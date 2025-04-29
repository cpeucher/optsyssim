% -------------------------------------------------------------------------
% Test of opt_tf_mrr_crow function
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
file_name_core_figure = strrep(mfilename,'test','fig');
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

line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));



% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
global CONSTANT
global reference_frequency

CONSTANT = core_load_constants();

reference_frequency = 193.1e12;
df = 10e6;
nsamples = 2^16;

sample_rate = nsamples*df;
dt = 1/sample_rate;
time_array = (0:nsamples-1)*dt;
frequency_array = (-nsamples/2:nsamples/2-1)*df;



%%
% -------------------------------------------------------------------------
% Single ring - all-pass configuration
% -------------------------------------------------------------------------

params_mrr.ng = 4.3;
% Waveguide group index
% Value for a silicon waveguide taken from
% W. Bogaerts et al., "Silicon microring resonators," Laser & Photonics 
% Reviews 6, 47 (2012) [DOI: 10.1002/lpor.201100017]



params_mrr.field_round_trip_loss = 0.85;
params_mrr.power_coupling_1 = 1 - 0.9^2;
params_mrr.power_coupling_2 = 0;
params_mrr.centre_frequency = 0;
params_mrr.fsr = 100e9;
vis_mrr.visualiser_status = 0;
vis_mrr.save_tf.status = 0;
vis_mrr.save_tf.file_name = 'mrr_tf.dat';
tf_sr = opt_tf_mrr(frequency_array,params_mrr,vis_mrr); 
% "Basic" transfer function of add-drop filter


ll = CONSTANT.c/params_mrr.ng/params_mrr.fsr;
% Ring circumference, in m
rr = ll/2/pi;
% Corresponding radius, assuming circular geometry, in m

params_mrr.loss_log = -10*log10(params_mrr.field_round_trip_loss^2)/ll;
% Ring waveguide loss, in dB/m


params_crow.nrings = 1;
params_crow.power_coupling = [params_mrr.power_coupling_1 0]';
params_crow.radius = rr;
params_crow.phase_shift = 0;
params_crow.neff = params_mrr.ng;
params_crow.loss_log = params_mrr.loss_log;
M = opt_tf_mrr_crow(frequency_array,params_crow);
% tf_through = squeeze(-M(1,1,:)./M(1,2,:)).';
% tf_drop = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';
tf_all_pass = squeeze((M(1,1,:) - M(2,1,:))./(M(2,2,:) - M(1,2,:))).';


fig_name = [file_name_core_figure '_single_all-pass'];
hfig = figure('Name',fig_name);
plot(frequency_array/1.0e9,10*log10(abs(tf_sr.s11).^2),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/1.0e9,10*log10(abs(tf_all_pass).^2),'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('relative frequency (GHz) ','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('analytical','scattering matrices','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])


% params_mrr = mrr_cmt_conversion(params_mrr);


%%
% -------------------------------------------------------------------------
% Single ring - all-pass configuration
% 
% We check the impact of loss and the critical coupling condition
% -------------------------------------------------------------------------


params_mrr.ng = 4.3;
% Waveguide group index
% Value for a silicon waveguide taken from
% W. Bogaerts et al., "Silicon microring resonators," Laser & Photonics 
% Reviews 6, 47 (2012) [DOI: 10.1002/lpor.201100017]


% Vary tt and check that the minimum transmission (at resonance) occurs
% when the round-trip field loss coefficient a is a = tt
tt = 0.85;
% tt = 0.9;
% tt = 0.98;
% tt = 0.999;

kappa = sqrt(1 - tt.^2);
% Cross coupling coefficient of the coupler (field)

 a = [0.6:0.01:1];
% Round-trip field loss coefficient

% a = [0.9:0.001:1];
% To compare between tt = 0.98 and tt = 0.999
% See Figure 2 in A. Yariv, "Universal relations for coupling of optical
% power between microresonators and dielectric waveguides," Electronics 
% Letters 36, 321 (2000) [DOI: 10.1049/el:20000340].



params_mrr.power_coupling_1 = kappa^2;
params_mrr.power_coupling_2 = 0;
params_mrr.centre_frequency = 0;
params_mrr.fsr = 100e9;
vis_mrr.visualiser_status = 0;
vis_mrr.save_tf.status = 0;
vis_mrr.save_tf.file_name = 'mrr_tf.dat';


ll = CONSTANT.c/params_mrr.ng/params_mrr.fsr;
% Ring circumference, in m
rr = ll/2/pi;
% Corresponding radius, assuming circular geometry, in m


params_crow.nrings = 1;
params_crow.power_coupling = [params_mrr.power_coupling_1 0]';
params_crow.radius = rr;
params_crow.phase_shift = 0;
params_crow.neff = params_mrr.ng;




tmin_sr = zeros(1,length(a));
tmin_crow = zeros(1,length(a));
% To store the minimum transmission at resonance

imin_sr = zeros(1,length(a));
imin_crow = zeros(1,length(a));
% We will check the corresponding indices


for ia = 1:length(a)

params_mrr.field_round_trip_loss = a(ia);

params_mrr.loss_log = -10*log10(params_mrr.field_round_trip_loss^2)/ll;
% Ring waveguide loss, in dB/m

tf_sr = opt_tf_mrr(frequency_array,params_mrr,vis_mrr); 
% "Basic" transfer function of add-drop filter

params_crow.loss_log = params_mrr.loss_log;
M = opt_tf_mrr_crow(frequency_array,params_crow);
% tf_through = squeeze(-M(1,1,:)./M(1,2,:)).';
% tf_drop = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';
tf_all_pass = squeeze((M(1,1,:) - M(2,1,:))./(M(2,2,:) - M(1,2,:))).';


tf_sr_1fsr = tf_sr.s11(abs(frequency_array) <= params_mrr.fsr/2);
tf_crow_1fsr = tf_all_pass(abs(frequency_array) <= params_mrr.fsr/2);
% For the sake of good order, we restrict the transfer function to one FSR.

[tmin_sr(ia), imin_sr(ia)] = min(abs(tf_sr.s11).^2);
[tmin_crow(ia), imin_crow(ia)] = min(abs(tf_all_pass).^2);

end








fig_name = [file_name_core_figure '_ap_critical-coupling'];
hfig = figure('Name',fig_name);
plot(a,tmin_sr,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(a,tmin_crow,'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
ylabel('transmission at resonance','Interpreter',fig.interpreter)
xlabel('round trip field loss $a$','Interpreter',fig.interpreter)
legend('analytical','numerical','Location','North','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])
vline(tt,'k--')



fig_name = [file_name_core_figure '_ap_critical-coupling_resonance_indices'];
figure('Name',fig_name)
plot(a,imin_sr,'b')
hold on
plot(a,imin_crow,'r--')
xlabel('a')
ylabel('resonance index')










%%
% -------------------------------------------------------------------------
% Single ring - add-drop configuration
% -------------------------------------------------------------------------

clear params_mrr
clear params_crow

params_mrr.ng = 4.3;
% Waveguide group index
% Value for a silicon waveguide taken from
% W. Bogaerts et al., "Silicon microring resonators," Laser & Photonics 
% Reviews 6, 47 (2012) [DOI: 10.1002/lpor.201100017]


params_mrr.field_round_trip_loss = 0.85;
params_mrr.power_coupling_1 = 1 - 0.9^2;
params_mrr.power_coupling_2 = params_mrr.power_coupling_1;
params_mrr.centre_frequency = 0;
params_mrr.fsr = 100e9;
vis_mrr.visualiser_status = 0;
vis_mrr.save_tf.status = 0;
vis_mrr.save_tf.file_name = 'mrr_tf.dat';
tf_sr = opt_tf_mrr(frequency_array,params_mrr,vis_mrr); 


ll = CONSTANT.c/params_mrr.ng/params_mrr.fsr;
% Ring circumference, in m
rr = ll/2/pi;
% Corresponding radius, assuming circular geometry, in m

params_mrr.loss_log = -10*log10(params_mrr.field_round_trip_loss^2)/ll;
% Ring waveguide loss, in dB/m


params_crow.nrings = 1;
params_crow.power_coupling = [params_mrr.power_coupling_1 params_mrr.power_coupling_2]';
params_crow.radius = rr;
params_crow.phase_shift = 0;
params_crow.neff = params_mrr.ng;
params_crow.loss_log = params_mrr.loss_log;
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through_2 = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop_2 = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';
% tf_all_pass = squeeze((M(1,1,:) - M(2,1,:))./(M(2,2,:) - M(1,2,:))).';


fig_name = [file_name_core_figure '_single_add-drop_through'];
hfig = figure('Name',fig_name);
plot(frequency_array/1.0e9,10*log10(abs(tf_sr.s11).^2),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/1.0e9,10*log10(abs(tf_through_2).^2),'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('relative frequency (GHz) ','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('analytical','scattering matrices','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])


fig_name = [file_name_core_figure '_single_add-drop_drop'];
hfig = figure('Name',fig_name);
plot(frequency_array/1.0e9,10*log10(abs(tf_sr.s21).^2),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/1.0e9,10*log10(abs(tf_drop_2).^2),'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('relative frequency (GHz) ','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('analytical','scattering matrices','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])


%%
% -------------------------------------------------------------------------
% Coupled rings: 2 rings 
%
% We reproduce some results from 
% T. Kato and Y. Kokubun, "Optimum coupling coefficients in second-order 
% series-coupled ring resonator for nonblocking wavelength channel switch," 
% Journal of Lightwave Technology 24, 991 (2006)
% [DOI: 10.1109/JLT.2005.862425].
%
% -------------------------------------------------------------------------
clear params_crow

% We can play with the parameters a, kappab, kappar

a = 0.9;
% Field round trip loss (linear)
kappab = 0.2;
% Power coupling coefficient between bus and ring
kappar = 0.012;
% Power coupling coefficient between rings

fsr = 100e9;
% Free spectral range, in Hz
% Value is arbitrary.

% For the matrix calculation, we need to introduce physical parameters
% including ring radii and indices
ng = 4.3;
% Arbitrary. c.f. supra.
ll = CONSTANT.c/ng/fsr;
% Ring circumference, in m
rr = ll/2/pi;
% Ring radius, in m

loss_log = -10*log10(a^2)/ll;
% Ring waveguide loss, in dB/m

params_crow.nrings = 2;
params_crow.power_coupling = [kappab kappar kappab]';
params_crow.radius = [rr rr];
params_crow.phase_shift = [0 0];
params_crow.neff = [ng ng]';
params_crow.loss_log = [loss_log loss_log]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';
% tf_all_pass = squeeze((M(1,1,:) - M(2,1,:))./(M(2,2,:) - M(1,2,:))).';


% We implement the analytical calculation according to eqs (7)-(8) in the
% reference by Kato and Kokubun
tb = 1 - kappab;
tr = 1 - kappar;


phi = 2*pi*frequency_array/fsr;
% Phase shift

z = exp(-1i*phi);
% Phase shift term

ht = (a*sqrt(tb)*z.^2 - sqrt(a)*(1 + tb)*sqrt(tr)*z + sqrt(tb))./(a*tb*z.^2 - 2*sqrt(a*tb*tr)*z + 1);
% In-to-through transfer function
% Eq. (7) in Kato and Kokubun 

hd = (sqrt(a)*(1 - tb)*sqrt(1 - tr)*z)./(a*tb*z.^2 - 2*sqrt(a*tb*tr)*z + 1);
% In-to-drop transfer function
% Eq. (8) in Kato and Kokubun

fig_name = [file_name_core_figure '_crow2_kato_kb' num2str(kappab)];
hfig = figure('Name',fig_name);
plot(frequency_array/fsr,10*log10(abs(tf_through).^2),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/fsr,10*log10(abs(tf_drop).^2),'Color','r','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(frequency_array/fsr,10*log10(abs(ht).^2),'Color','k','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(frequency_array/fsr,10*log10(abs(hd).^2),'Color','k','LineStyle','-.','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('$\left(\nu - \nu_c\right)/\Delta\nu_\mathrm{FSR}$','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('through - numerical','drop - numerical','through - analytical','drop - analytical','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([-0.2 0.2])
ylim([-60 15])
text(-0.18,-3,sprintf('$$a$$ = %3.2f',a),'Interpreter',fig.interpreter,'fontsize',18);
text(-0.18,-10,sprintf('$$\\kappa_b$$ = %3.2f',kappab),'Interpreter',fig.interpreter,'fontsize',18);
text(-0.18,-17,sprintf('$$\\kappa_r$$ = %3.4f',kappar),'Interpreter',fig.interpreter,'fontsize',18);


% Comparison of transfer functions for kb = 0.12, 0.012 and 0.0012

a = 1;
loss_log = -10*log10(a^2)/ll;
params_crow.loss_log = [loss_log loss_log]';

kappab = 0.2;



kappar1 = 0.12;
tr1 = 1 - kappar1;
params_crow.power_coupling = [kappab kappar1 kappab]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through1 = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop1 = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';

ht1 = (a*sqrt(tb)*z.^2 - sqrt(a)*(1 + tb)*sqrt(tr1)*z + sqrt(tb))./(a*tb*z.^2 - 2*sqrt(a*tb*tr1)*z + 1);
hd1 = (sqrt(a)*(1 - tb)*sqrt(1 - tr1)*z)./(a*tb*z.^2 - 2*sqrt(a*tb*tr1)*z + 1);


kappar2 = 0.012;
tr2 = 1 - kappar2;
params_crow.power_coupling = [kappab kappar2 kappab]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through2 = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop2 = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';

ht2 = (a*sqrt(tb)*z.^2 - sqrt(a)*(1 + tb)*sqrt(tr2)*z + sqrt(tb))./(a*tb*z.^2 - 2*sqrt(a*tb*tr2)*z + 1);
hd2 = (sqrt(a)*(1 - tb)*sqrt(1 - tr2)*z)./(a*tb*z.^2 - 2*sqrt(a*tb*tr2)*z + 1);


kappar3 = 0.0012;
tr3 = 1 - kappar3;
params_crow.power_coupling = [kappab kappar3 kappab]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through3 = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop3 = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';

ht3 = (a*sqrt(tb)*z.^2 - sqrt(a)*(1 + tb)*sqrt(tr3)*z + sqrt(tb))./(a*tb*z.^2 - 2*sqrt(a*tb*tr3)*z + 1);
hd3 = (sqrt(a)*(1 - tb)*sqrt(1 - tr3)*z)./(a*tb*z.^2 - 2*sqrt(a*tb*tr3)*z + 1);


fig_name = [file_name_core_figure '_crow2_kato_kb' num2str(kappab) '_through'];
hfig = figure('Name',fig_name);
plot(frequency_array/fsr,10*log10(abs(tf_through1).^2),'Color',line_color(1,:),'LineStyle','-.','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/fsr,10*log10(abs(tf_through2).^2),'Color',line_color(2,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(frequency_array/fsr,10*log10(abs(tf_through3).^2),'Color',line_color(3,:),'LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('$\left(\nu - \nu_c\right)/\Delta\nu_\mathrm{FSR}$','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('$\kappa_r = 0.12$','$\kappa_r = 0.012$','$\kappa_r = 0.0012$','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([-0.2 0.2])
ylim([-35 15])
text(-0.18,-3,sprintf('$$a$$ = %3.2f',a),'Interpreter',fig.interpreter,'fontsize',18);
text(-0.18,-10,sprintf('$$\\kappa_b$$ = %3.2f',kappab),'Interpreter',fig.interpreter,'fontsize',18);
text(-0.18,-17,sprintf('$$\\kappa_r$$ = %3.4f',kappar),'Interpreter',fig.interpreter,'fontsize',18);


fig_name = [file_name_core_figure '_crow2_kato_kb' num2str(kappab) '_drop'];
hfig = figure('Name',fig_name);
plot(frequency_array/fsr,10*log10(abs(tf_drop1).^2),'Color',line_color(1,:),'LineStyle','-.','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/fsr,10*log10(abs(tf_drop2).^2),'Color',line_color(2,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(frequency_array/fsr,10*log10(abs(tf_drop3).^2),'Color',line_color(3,:),'LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('$\left(\nu - \nu_c\right)/\Delta\nu_\mathrm{FSR}$','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('$\kappa_r = 0.12$','$\kappa_r = 0.012$','$\kappa_r = 0.0012$','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xlim([-0.5 0.5])
ylim([-60 15])
text(-0.45,-3,sprintf('$$a$$ = %3.2f',a),'Interpreter',fig.interpreter,'fontsize',18);
text(-0.45,-10,sprintf('$$\\kappa_b$$ = %3.2f',kappab),'Interpreter',fig.interpreter,'fontsize',18);
text(-0.45,-17,sprintf('$$\\kappa_r$$ = %3.4f',kappar),'Interpreter',fig.interpreter,'fontsize',18);

% Conclusion
% We observe a good match between the analytical transfer functions
% established in Kato & Kokubun and the numerical results obtained
% through the numerical scattering matrices approach.
% We reproduce Figures 3 and 5 in Kato & Kokubun with our numerical 
% scattering matrices approach



%%
% -------------------------------------------------------------------------
% Coupled rings: 6 rings
% We reproduce the results (Figure 4) of
% R. Orta, P. Savi, R. Tascone, and D. Trinchero, "Synthesis of multiple-
% ring-resonator filters for optical systems," IEEE Photonics Technology 
% Letters 7, 1447 (1995) [DOI: 10.1109/68.477278].
%
% This example is also treated in:
% C. K. Madsen and J. H. Zhao, Optical filter design and analysis-  a
% signal processing approach, Wiley, 1999
% See Figure 5.15, page 258.
% -------------------------------------------------------------------------


fsr = 100e9;
% Free spectral range, in Hz
% Value is arbitrary.

% For the matrix calculation, we need to introduce physical parameters
% including ring radii and indices
ng = 4.3;
% Arbitrary. c.f. supra.
ll = CONSTANT.c/ng/fsr;
% Ring circumference, in m
rr = ll/2/pi;
% Ring radius, in m


% The couplers are specified in terms of coupling angle in Orta et al.
phik = [0.2064 0.0524 0.0347 0.0323 0.0347 0.0524 0.2064]*pi;
% In order to obtain the cross power coupling coefficients listed in Fig.
% 5.15 of Madsen and Zhao, we need to take:
params_crow.power_coupling = (sin(phik).^2)';

params_crow.nrings = 6;
params_crow.radius = rr;
params_crow.phase_shift = 0;
params_crow.neff = ng;
params_crow.loss_log = 0;
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';


fig_name = [file_name_core_figure '_madsen_fig5-15'];
hfig = figure('Name',fig_name);
plot(frequency_array/fsr,10*log10(abs(tf_through).^2),'Color',line_color(1,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/fsr,10*log10(abs(tf_drop).^2),'Color',line_color(2,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('$\left(\nu - \nu_c\right)/\Delta\nu_\mathrm{FSR}$','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('through','drop','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter; 
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.YTick = [-140:20:0];
ax.XTick = [-0.5:0.2:0.5];
ax.XMinorTick = 'on';
xlim([-0.5 0.5])
ylim([-140 5])

% We indeed reproduce Fig. 5.15 in Madsen & Zhao. 
% The 20-dB bandwidth is 0.082xFSR, as specified in the text.
% The figure, however, does not match Figure 4 of Orta et al.




%%
% -------------------------------------------------------------------------
% Coupled rings: 6 rings
% We attempt to reproduce the results of Figure 5.16
% C. K. Madsen and J. H. Zhao, Optical filter design and analysis-  a
% signal processing approach, Wiley, 1999
% We observe we had a good match with Fig. 5.15.
%
% -------------------------------------------------------------------------

fsr = 100e9;
% Free spectral range, in Hz
% Value is arbitrary.

% For the matrix calculation, we need to introduce physical parameters
% including ring radii and indices
ng = 4.3;
% Arbitrary. c.f. supra.
ll = CONSTANT.c/ng/fsr;
% Ring circumference, in m
rr = ll/2/pi;
% Ring radius, in m

params_crow.power_coupling = [0.8465 0.2753 0.0940 0.0695 0.0940 0.2753 0.8465]';
% Coupling parameters as specified in Fig. 5-16 of Madsen and Zhao

params_crow.nrings = 6;
params_crow.radius = rr;
params_crow.phase_shift = 0;
params_crow.neff = ng;
params_crow.loss_log = 0;
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';


fig_name = [file_name_core_figure '_madsen_fig5-16'];
hfig = figure('Name',fig_name);
plot(frequency_array/fsr,10*log10(abs(tf_through).^2),'Color',line_color(1,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(frequency_array/fsr,10*log10(abs(tf_drop).^2),'Color',line_color(2,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('$\left(\nu - \nu_c\right)/\Delta\nu_\mathrm{FSR}$','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('through','drop','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter; 
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.YTick = [-140:20:0];
ax.XTick = [-0.5:0.2:0.5];
ax.XMinorTick = 'on';
xlim([-0.5 0.5])
ylim([-140 5])

% The 20-dB bandwidth is 0.24xFSR, as specified in the text.





%% 
% -------------------------------------------------------------------------
% Coupled rings
% 
% B.E. Little, et al., "Microring resonator channel dropping filters," 
% Journal of Lightwave Technology 15, 998 (1997) [DOI: 10.1109/50.588673].
% 
% -------------------------------------------------------------------------

% clear all
% CONSTANT = core_load_constants();

reference_frequency = CONSTANT.c/1334e-9;
% We change the reference frequency, to match the paper.

% We also need a large simulation bandwidth: 
df = 20e6;
nsamples = 2^16;
sample_rate = nsamples*df;
dt = 1/sample_rate;
time_array = (0:nsamples-1)*dt;
frequency_array = (-nsamples/2:nsamples/2-1)*df;


lambda = CONSTANT.c./(reference_frequency + frequency_array);
% Wavelength axis, in nm


R = 1.7e-6;
ng = 3;
% Ring parameters.

kappa = 0.3;
% Dunno what it should be. We just quickly tune it to approximately match 
% Fig. 5(b) in the paper.

params_crow.radius = R;
params_crow.phase_shift = 0;
params_crow.neff = ng;
params_crow.loss_log = 0;
% This parameters are identical, regardless of the number of rings.


params_crow.nrings = 1;
params_crow.power_coupling = [kappa^2 kappa^2]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through_1 = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop_1 = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';


kappa1 = sqrt(0.250)*kappa^2;
% Coupling coefficient for N = 2 rings
% Table I in the paper

params_crow.nrings = 2;
params_crow.power_coupling = [kappa^2 kappa1^2 kappa^2]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through_2 = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop_2 = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';



kappa1 = sqrt(0.125)*kappa^2;
% Coupling coefficient for N = 3 rings
% Table I in the paper



params_crow.nrings = 3;
params_crow.power_coupling = [kappa^2 kappa1^2 kappa1^2 kappa^2]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through_3 = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop_3 = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';



fig_name = [file_name_core_figure '_crow_little_drop'];
hfig = figure('Name',fig_name);
plot(lambda/1.0e-9,10*log10(abs(tf_drop_1).^2),'Color',line_color(1,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(lambda/1.0e-9,10*log10(abs(tf_drop_2).^2),'Color',line_color(2,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(lambda/1.0e-9,10*log10(abs(tf_drop_3).^2),'Color',line_color(3,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')

xlabel('wavelength (nm) ','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('1 ring','2 rings','3 rings','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([1330 1340])
ylim([-80 0])



fig_name = [file_name_core_figure '_crow_little_through'];
hfig = figure('Name',fig_name);
plot(lambda/1.0e-9,10*log10(abs(tf_through_1).^2),'Color',line_color(1,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(lambda/1.0e-9,10*log10(abs(tf_through_2).^2),'Color',line_color(2,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(lambda/1.0e-9,10*log10(abs(tf_through_3).^2),'Color',line_color(3,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')

xlabel('wavelength (nm) ','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('1 ring','2 rings','3 rings','Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([1330 1340])
ylim([-80 0])




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

