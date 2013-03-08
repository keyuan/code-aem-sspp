% Code for reproducing figure 2 in the manuscript
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
option = inferset('model','sspp','method','em','maxiter',500,   ...
    'tolfun',1e-6,'fixparam',{'sigma'},'intype','spike','stadim', ...
    1,'estep','approxsmoother','display','iter','fltopt','newton'...
    ,'cif','exp');

%% Set dimensions 
dim         = struct('type','dimension of data' );
dim.delta   = 1e-2;                             % Time resolution                               
dim.totchan = 10;                               % Number of channels
dim.tottime = 20;                               % Total observation time
dim.stadim  = 1;                                % Dimension of state
dim.inparam = [1,ceil(1/dim.delta)];            % Input parameter

%% Set parameters 
param = struct('type','parameters of sspp');
param.true.rho     = 0.6;                       % AR coeffeicient
param.true.alpha   = 4;                         % Inpute weight
param.true.sigmasq = 0.01;                      % State noise variance                  
param.true.mu      = 0;                         % Background firing rate
param.true.beta    = ones(1,dim.totchan);       % State weight
param.true.xinit   = 0;                         % Initial state
param.true.covinit = param.true.sigmasq/...     % Initial variance
                     (1-param.true.rho^2);        
param.true.gamma   = [];                        % History coefficient

%% Generate synthetic data 
use_orginal_data = 1;
if (~use_orginal_data)
    synthdata = synthdatapp(dim, param, option);
end
%% Save data
if (~use_orginal_data)
    savedata = 0;
    if (savedata)        
        datafilepath = './data/synthdb/';
        filename     = ['data_fig_2_new' '.mat']; 
        save([datafilepath filename],'synthdata')
    end
end
%% Inference and learning via EM 
% load data
if (use_orginal_data)
    load('./data/synthdb/data_fig_2.mat')
end

%% Set initial conditions
param.est.rho     = 0.1;
param.est.alpha   = 4;
param.est.sigmasq = param.true.sigmasq;
param.est.mu      = 1;
param.est.beta    = 2*param.true.beta;
param.est.xinit   = rand;
param.est.covinit = 1;

%% EM
tic
[param,stats,lbsave,nem] = em_sspp(synthdata,param,option);
runtime = toc; 

%% Draw figures 
showresult = 1;
if (showresult)

set(0,'defaulttextinterpreter','latex');
red = [0.9, 0, 0];
blue = [0, 0, 0.7];
green = [0, 0.7, 0];


figure(2),clf 
subplot(232)
plot(param.save.rho(1:nem-1),'color',blue,'linewidth',1.5)
hold on 
plot(param.true.rho(:,ones(1,nem-1)),'linestyle','--','color', blue, ...
    'linewidth',1.5)
hold off
xlim([0,nem])
ylim([0,1])
xlabel('Iteration')   
ylabel('$\rho$')


subplot(233)
plot(param.save.alpha(:,1:nem-1)','linewidth',1.5)
hold on 
plot(param.true.alpha(ones(1,nem-1),:),'--','linewidth',1.5)
hold off
xlim([0,nem])
ylim([2,4.5])
xlabel('Iteration')   
ylabel('$\alpha$')

subplot(235)
plot(param.save.mu(1:nem-1),'color',blue,'linewidth',1.5)
hold on 
plot(param.true.mu(:,ones(1,nem-1)),'linestyle','--','color', blue, ...
    'linewidth',1.5)
hold off
xlim([0,nem])
ylim([-0.2,1])
xlabel('Iteration')   
ylabel('$\mu$')  

subplot(236)
plot(param.save.beta(:,1:nem-1)','linewidth',1.5)
hold on 
plot(param.true.beta(:,ones(1,nem-1))','--','linewidth',1.5)
hold off
xlim([0,nem])
ylim([0.8,2])
xlabel('Iteration')   
ylabel('$\beta_c$') 

subplot(2,3,[1,4])
plot(lbsave(2:nem),'color',blue,'linewidth',1.5)
xlim([0,nem])
xlabel('Iteration')   
ylabel('$\mathcal{Q}(\theta,q)$')

% % Creat a zoomin subfigure
% figure(1),clf
% plot(lbsave(2:nem),'color',blue,'linewidth',1.5)
% xlim([0,nem])
% ylim([290,330])



set(gcf,'units','centimeters');
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2),15,10]);

matlabfrag('./fig/fig-2')
end



