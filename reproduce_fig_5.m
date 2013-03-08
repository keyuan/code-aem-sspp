% Reproduce figure 5 in the manuscript
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
option = inferset('model','sspp','method','em','maxiter',50,    ...
    'tolfun',1e-6,                                               ...
    'intype','spike','stadim',                                   ...
    1,'estep','approxsmoother','display','iter','fltopt','newton'...
    ,'cif','exp');

%% Set dimensions 
dim         = struct('type','dimension of data' );
dim.delta   = 1e-2;                             % Time resolution                               
dim.totchan = 1;                                % Number of channels
dim.tottime = 20;                               % Total observation time
dim.stadim  = 1;                                % Dimension of state
dim.inparam = [1,ceil(1/dim.delta)];            % Input parameter

%% Set parameters 
param = struct('type','parameters of sspp');
param.true.rho     = 0.8;                       % AR coeffeicient
param.true.alpha   = 4;                         % Inpute weight
param.true.sigmasq = 0.01;                      % State noise variance                  
param.true.mu      = 0;                         % Background firing rate
param.true.beta    = ones(1,dim.totchan);       % State weight
param.true.xinit   = 0;                         % Initial state
param.true.covinit = param.true.sigmasq/...     % Initial variance
                     (1-param.true.rho^2);        
param.true.gamma   = [];                        % History coefficient

%% Generate synthetic data 
synthdata = [];
% synthdata = synthdatapp(dim, param, option);

%% Inference and learning via EM 
% load data
if isempty(synthdata)
    load('./data/synthdb/data_synthdb_c_1_alpha_4.mat')
end

param.est.rho     = param.true.rho;
param.est.alpha   = param.true.alpha;
param.est.sigmasq = param.true.sigmasq;
param.est.mu      = param.true.mu;
param.est.beta    = param.true.beta(:);
param.est.xinit   = param.true.xinit;
param.est.covinit = param.true.covinit;

%% Compute Q-function 
% Get options ------------------------------------------------------------%
estep    = inferget(option,'estep');
tol      = inferget(option,'tolfun');
intype   = inferget(option,'intype');
stadim   = inferget(option,'stadim');
totem    = inferget(option,'maxiter');
display  = inferget(option,'display');
%-------------------------------------------------------------------------%

% Get data ---------------------------------------------------------------%
y      = synthdata.y;
time   = synthdata.time;
delta  = synthdata.delta;
if ~isempty(synthdata.in)
    in = synthdata.in;
else
    inparam = synthdata.inparam;
    in   = input_sspp(time,intype,inparam);
end
%-------------------------------------------------------------------------%

% Get dimensions ---------------------------------------------------------%
totsamp = size(y,2);
lbsave  = zeros(1,totem);
%-------------------------------------------------------------------------%

% Preallocate saved parameter estimaiton result --------------------------%
if stadim == 1;  
    param.save.rho    = zeros(stadim,totem);
    param.save.alpha  = zeros(length(param.est.alpha),totem);
    param.save.mu     = zeros(1,totem);
    param.save.beta   = zeros(size(y,1),totem);
else
    param.save.rho    = zeros(stadim,stadim,totem);
    param.save.alpha  = zeros(stadim,length(param.est.alpha),totem);
    param.save.mu     = zeros(size(y,1),totem);
    param.save.beta   = zeros(stadim,size(y,1),totem);
end
%-------------------------------------------------------------------------%

% Align initial guesses --------------------------------------------------%
if stadim == 1    
    param.save.rho(:,1)     = param.est.rho;
    param.save.alpha(:,1)   = param.est.alpha;
    param.save.mu(:,1)      = param.est.mu;
    param.save.beta(:,1)    = param.est.beta(:);
else 
    param.save.rho(:,:,1)   = param.est.rho;
    param.save.alpha(:,:,1) = param.est.alpha;
    param.save.mu(:,1)      = param.est.mu;
    param.save.beta(:,:,1)  = param.est.beta;
end
%-------------------------------------------------------------------------%

% Preallocate E-step variables -------------------------------------------%
stats          = struct('type','inferred statistics');
if stadim == 1
    stats.xpred    = zeros(stadim,totsamp);
    stats.covpred  = zeros(stadim,totsamp);
    stats.xpost    = zeros(stadim,totsamp);
    stats.covpost  = zeros(stadim,totsamp);
    stats.xsmth    = zeros(stadim,totsamp);
    stats.covsmth  = zeros(stadim,totsamp);
    stats.autoexp  = zeros(stadim,totsamp);
    stats.crosscov = zeros(stadim,totsamp-1);
    stats.crossexp = zeros(stadim,totsamp-1);
else
    stats.xpred    = zeros(stadim,totsamp);
    stats.covpred  = zeros(stadim,stadim,totsamp);
    stats.xpost    = zeros(stadim,totsamp);
    stats.covpost  = zeros(stadim,stadim,totsamp);
    stats.xsmth    = zeros(stadim,totsamp);
    stats.covsmth  = zeros(stadim,stadim,totsamp);
    stats.autoexp  = zeros(stadim,stadim,totsamp);
    stats.crosscov = zeros(stadim,stadim,totsamp-1);
    stats.crossexp = zeros(stadim,stadim,totsamp-1);
end
%-------------------------------------------------------------------------%

%% qfunc 1
num = 100;
rho1 = linspace(-1,1,num);
alpha1 = linspace(0,6,num);
sigmasq1 = linspace(0,0.02,num);
mu1 = linspace(-5,5,num);
beta1 = linspace(-10,5,num)';
beta1 = repmat(beta1,1,dim.totchan);
qfunc1 = zeros(5,num);

for jj = 1:5
    param.est.rho     = param.true.rho;
    param.est.alpha   = param.true.alpha;
    param.est.sigmasq = param.true.sigmasq;
    param.est.mu      = param.true.mu;
    param.est.beta    = param.true.beta;
    param.est.xinit   = param.true.xinit;
    param.est.covinit = param.true.covinit;
    for ii = 1:num
        [jj,ii]
        switch jj
            case 1, param.est.rho = rho1(ii);
            case 2, param.est.alpha = alpha1(ii);
            case 3, param.est.sigmasq = sigmasq1(ii);
            case 4, param.est.mu = mu1(ii);
            case 5, param.est.beta = beta1(ii,:);
        end
        [none,lb] = feval(estep, y, in, delta, stats, param, option); 
        qfunc1(jj,ii) = lb;
    end
end

%%

red = [0.9, 0, 0];
blue = [0, 0, 0.7];
green = [0, 0.7, 0];

true(1) = param.true.rho;
true(2) = param.true.alpha;
true(3) = param.true.sigmasq;
true(4) = param.true.mu;
true(5) = mean(param.true.beta);

figure(1)
set(0,'defaulttextinterpreter','latex');
scale = 100;
subplot(521)
plot(rho1,qfunc1(1,:)/scale)
[ma,ind] = max(qfunc1(1,:)); mi = min(qfunc1(1,:));
hold on 
plot([rho1(ind); rho1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(1); true(1)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([15 16.34])
xlabel('$\rho$')
ylabel(['$\mathcal{Q}(q,\rho)$' ...
    ' ($\times ' num2str(scale) '$)'])
hold off

subplot(522)
plot(rho1,qfunc1(1,:)/scale)
[ma,ind] = max(qfunc1(1,:)); mi = min(qfunc1(1,:));
hold on 
plot([rho1(ind); rho1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(1); true(1)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([16.2 16.34])
xlabel('$\rho$')
hold off


subplot(523)
plot(alpha1,qfunc1(2,:)/scale)
hold on 
[ma,ind] = max(qfunc1(2,:)); mi = min(qfunc1(2,:));
plot([alpha1(ind); alpha1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(2); true(2)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([15.8 16.34])
xlabel('$\alpha$')
ylabel(['$\mathcal{Q}(q,\alpha)$' ...
    ' ($\times ' num2str(scale) '$)'])
hold off

subplot(524)
plot(alpha1,qfunc1(2,:)/scale)
hold on 
[ma,ind] = max(qfunc1(2,:)); mi = min(qfunc1(2,:));
plot([alpha1(ind); alpha1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(2); true(2)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([16.2 16.34])
xlim([2.5,5])
xlabel('$\alpha$')
hold off

subplot(525)
plot(sigmasq1,qfunc1(3,:)/scale)
hold on 
[ma,ind] = max(qfunc1(3,:)); mi = min(qfunc1(3,:));
plot([sigmasq1(ind); sigmasq1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(3); true(3)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([12,49])
xlabel('$\sigma^2_{\varepsilon}$')
ylabel(['$\mathcal{Q}(q,\sigma^2_{\varepsilon})$' ...
    ' ($\times ' num2str(scale) '$)'])
hold off

subplot(526)
plot(sigmasq1,qfunc1(3,:)/scale)
hold on 
[ma,ind] = max(qfunc1(3,:)); mi = min(qfunc1(3,:));
plot([sigmasq1(ind); sigmasq1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(3); true(3)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([20,50])
xlabel('$\sigma^2_{\varepsilon}$')
hold off

subplot(527)
plot(mu1,qfunc1(4,:)/scale)
hold on
[ma,ind] = max(qfunc1(4,:)); mi = min(qfunc1(4,:));
plot([mu1(ind); mu1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(4); true(4)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([15 16.5])
xlim([-5,2])
xlabel('$\mu$')
ylabel(['$\mathcal{Q}(q,\mu)$' ...
    ' ($\times ' num2str(scale) '$)'])
hold off

subplot(528)
plot(mu1,qfunc1(4,:)/scale)
hold on
[ma,ind] = max(qfunc1(4,:)); mi = min(qfunc1(4,:));
plot([mu1(ind); mu1(ind)],[mi;ma]/scale,'color',blue)
hold on
plot([true(4); true(4)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([16.2 16.35])
xlim([-1,0.6])
xlabel('$\mu$')
hold off

subplot(529)
plot(beta1(:,1)',qfunc1(5,:)/scale)
hold on
[ma,ind] = max(qfunc1(5,:)); mi = min(qfunc1(5,:));
plot([beta1(ind,1); beta1(ind,1)],[mi;ma]/scale,'color',blue)
hold on
plot([true(5); true(5)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([15 16.5])
xlim([-1,1.8])
xlabel('$\beta_c$')
ylabel(['$\mathcal{Q}(q,\beta_c)$' ...
    ' ($\times ' num2str(scale) '$)'])
hold off

subplot(5,2,10)
plot(beta1(:,1)',qfunc1(5,:)/scale)
hold on
[ma,ind] = max(qfunc1(5,:)); mi = min(qfunc1(5,:));
plot([beta1(ind,1); beta1(ind,1)],[mi;ma]/scale,'color',blue)
hold on
plot([true(5); true(5)],[mi;ma]/scale,'color',green, 'linestyle', '--')
ylim([16.2 16.35])
xlim([0.6,1.3])
xlabel('$\beta_c$')
hold off

set(gcf,'units','centimeters');
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2),15,20]);

matlabfrag('./fig/fig-5')
