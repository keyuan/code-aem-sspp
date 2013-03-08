% Reproduce figure 1b in the manuscript
%
% -------------------------------------
% author: ke yuan 
% email : ky08r@ecs.soton.ac.uk

% system settings
clear all; close all; clc; 
randn('state',sum(100*clock));
rand('twister',sum(100*clock));
echo off;

%% Set path
pathset;

%% Set option
option = inferset('model','sspp','cif','exp','intype','spike');

%% Set dimensions 
dim         = struct('type','dimension of data');
dim.delta   = 1e-2;                                % Time resolution                               
dim.totchan = 10;                                  % Number of channels
dim.tottime = 20;                                  % Total observation time
dim.stadim  = 1;                                   % Dimension of state
dim.inparam = [1 ceil(1/dim.delta)];               % Input parameter

%% Set parameters 
param = struct('type','parameters of sspp');
param.true.rho     = 0.8;                          % AR coeffeicient
param.true.alpha   = 4;                            % Inpute weight
param.true.sigmasq = 0.01;                         % State noise variance                  
param.true.mu      = 0;                            % Background firing rate
param.true.beta    = linspace(0.5,1,...
                     dim.totchan);                 % State weight
param.true.xinit   = 0;                            % Initial state
param.true.covinit = param.true.sigmasq/...        % Initial variance
                     (1-param.true.rho^2);  
param.true.gamma   = [];
    
%% Generate synthetic data 
synthdata = synthdatapp(dim, param, option);

%% Draw figure
time  = synthdata.time;
y     = synthdata.y;
in    = synthdata.in;
xtrue = synthdata.x;

ciftype = inferget(option,'cif','default','fast');
if strcmp(ciftype,'sigmoid')
    cifd = synthdata.cif;
else cifd = synthdata.cif*dim.delta;
end

set(0,'defaulttextinterpreter','latex');

red = [0.9, 0, 0];
blue = [0, 0, 0.7];
green = [0, 0.7, 0];

figure(1), clf

idx = find(in==1);
u = 0.5;
lower = 5*ones(size(idx));
upper = 5+u*ones(size(idx));
plot([time(idx);time(idx)],[lower;upper],'color','k')
hold on 

plot(time,xtrue,'color',blue)
hold on

for ii = 1:size(y,1)
    idx = find(y(ii,:)==1);
    u = 0.5;
    lower = -1-(ii)*u*ones(size(idx));
    upper = -1-(ii+0.8)*u*ones(size(idx));
    plot([time(idx);time(idx)],[lower;upper],'color',green)
    hold on
end
ylim([-7,6])
xlabel('Time (s)')
ylabel('Data')

hold off

set(gcf,'units','centimeters');
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2),8,8]);

matlabfrag('./fig/fig-1b')


