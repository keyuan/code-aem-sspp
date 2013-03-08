% Reproduce figure 7 in the manuscript 
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
option = inferset('model','sspp','method','em','maxiter',600,   ...
    'tolfun',1e-6,'intype','spike','stadim',2, ...
    'estep','approxsmoother','display','iter','fltopt','fixpt'...
    ,'cif','exp');

%% Set dimensions 
dim         = struct('type','dimension of data' );
dim.delta   = 1e-2;                             % Time resolution                               
dim.totchan = 10;                               % Number of channels
dim.tottime = 20;                               % Total observation time
dim.stadim  = 2;                                % Dimension of state
dim.inparam = [3,ceil(1/dim.delta),...
    ceil(1.5/dim.delta),ceil(2/dim.delta)];     % Input parameter

%% Set parameters 
param = struct('type','parameters of sspp');
param.true.rho     = [0.8 0; -0.2 0.9];         % AR coeffeicient
param.true.alpha   = [.5 2 1.2; 1.1 1.3 1.19];  % Inpute weight
param.true.sigmasq = diag([0.01,0.01]);         % State noise variance                  
param.true.mu      = 0;                         % Background firing rate
beta               = linspace(0.5,1,dim.totchan);
beta               = repmat(beta,2,1);
param.true.beta    = beta;                      % State weight
param.true.xinit   = zeros(dim.stadim,1);       % Initial state
param.true.covinit = eye(dim.stadim);       
param.true.gamma   = [];                        % History coefficient

%% Generate synthetic data 
use_orginal_data = 1;
if (~use_orginal_data)
    synthdata = synthdatapp(dim, param, option);
end

%% Save data
savedata = 0;
if (savedata)        
    datafilepath = './data/synthdb/';
    filename     = ['data_fig_7_new' '.mat']; 
    save([datafilepath filename],'synthdata')
end

%% Inference and learning via EM 
% load data
if (use_orginal_data)
    load('./data/synthdb/data_fig_7.mat')
end

%% Fixed sigma and beta
option = inferset(option,'fixparam',{'sigma','beta'});

%% Set initial conditions
param.est.rho     = diag([0.2, 0.2]);
param.est.alpha   = rand(2,3);
param.est.sigmasq = param.true.sigmasq;
param.est.mu      = rand;
param.est.beta    = param.true.beta;
param.est.xinit   = rand(dim.stadim,1);
param.est.covinit = rand*eye(dim.stadim);

%% EM
tic
[param,stats,lbsave,nem] = em_sspp(synthdata,param,option);
runtime = toc; 

%% Draw figures 
showresult = 0;


set(0,'defaulttextinterpreter','latex');
red = [0.9, 0, 0];
blue = [0, 0, 0.7];
green = [0, 0.7, 0];

figure(1),clf

time  = synthdata.time;
xtrue = synthdata.x;
y     = synthdata.y;
in    = synthdata.in;
fixparam = inferget(option,'fixparam');
stadim   = inferget(option,'stadim');
totsamp = size(xtrue,2);

cov = zeros(stadim,totsamp);
for nsamp = 1:totsamp
    cov(:,nsamp) = diag(stats.covsmth(:,:,nsamp));
end

for ii = 1:stadim
    subplot(4,2,[ii,ii+2])
    plot(time,stats.xsmth(ii,:),'-','color',red)
    hold on
    plot(time,stats.xsmth(ii,:)+1.96*sqrt(cov(ii,:)),'--','color',red)
    hold on
    plot(time,stats.xsmth(ii,:)-1.96*sqrt(cov(ii,:)),'--','color',red)
    hold on 
    plot(time,xtrue(ii,:),'color',blue)
    hold on
    
    indim = size(in,1);
    for mm = 1:indim
        idx = find(in(mm,:)==1);
        u = 0.5;base = 4;
        lower = base+(mm)*u*ones(size(idx));
        upper = base+(mm+0.8)*u*ones(size(idx));
        plot([time(idx);time(idx)],[lower;upper],'color','k')
        hold on
    end

    for cc = 1:size(y,1)
        idx = find(y(cc,:)==1);
        u = 0.5;
        lower = -1-(cc)*u*ones(size(idx));
        upper = -1-(cc+0.8)*u*ones(size(idx));
        plot([time(idx);time(idx)],[lower;upper],'color',green)
        hold on
    end
    ylim([-7,7])
    xlabel('Time (s)')
    ylabel('Data and states esitmates')
    hold off
end

subplot(425)
plot(lbsave(:,2:nem-1)','color',blue,'linewidth',1.5)
xlabel('Iteration')   
ylabel('$\mathcal{Q}(q,\theta)$')


subplot(426)
plot(param.save.rho(:,1:nem-1)','linewidth',1.5)
hold on 
plot(repmat(param.true.rho(:),1,nem-1)','linestyle','--','linewidth',1.5)
hold off
xlabel('Iteration')   
ylabel('$\rho$')

subplot(427)
plot(param.save.alpha(:,1:nem-1)','linewidth',1.5)
hold on 
plot(repmat(param.true.alpha(:),1,nem-1)','--','linewidth',1.5)
hold off
ylim([0,2.1])
xlabel('Iteration')   
ylabel('$\alpha$')

subplot(428)
plot(param.save.mu(:,1:nem-1),'color',blue,'linewidth',1.5)
hold on 
plot(repmat(param.true.mu(:),1,nem-1),'color',blue,'linestyle','--', ...
    'linewidth',1.5)
hold off
ylim([-0.5,0.5])
xlabel('Iteration')   
ylabel('$\mu$')  

set(gcf,'units','centimeters');
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2),15,20]);

matlabfrag('./fig/fig-7')


