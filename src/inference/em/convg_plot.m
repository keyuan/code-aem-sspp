function [fig1] = convg_plot(data,stats,param,nem,lbsave,option)

time  = data.time;
xtrue = data.x;
y     = data.y;
in    = data.in;
fixparam = inferget(option,'fixparam');


set(0,'defaulttextinterpreter','latex');

red = [0.9, 0, 0];
blue = [0, 0, 0.7];
green = [0, 0.7, 0];

fig1 = figure;
subplot(4,2,1:4)
%= statss estimaiton results
idx = find(in==1);
u = 0.5;
lower = 5*ones(size(idx));
upper = 5+u*ones(size(idx));
plot([time(idx);time(idx)],[lower;upper],'color','k')
hold on


plot(time,stats.xsmth,'-','color',red)
hold on
plot(time,stats.xsmth+1.96*sqrt(stats.covsmth),...
     '--','color',red)
hold on
plot(time,stats.xsmth-1.96*sqrt(stats.covsmth),...
     '--','color',red)
hold on 

plot(time,xtrue,'-','color',blue)
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
ylabel('Data and states estimates')

hold off

% fig2 = figure(2);
%= free energy convergence performance
subplot(4,2,5)
plot(lbsave(2:nem),'color',blue,'linewidth',1.5)
ylim([-10000 1000])
xlim([0,nem])
xlabel('Iteration')   
ylabel('$\mathcal{Q}(\theta,q)$')

%= parameter convergence performance
subplot(4,2,6)
plot(param.save.rho(1:nem-1),'color',blue,'linewidth',1.5)
hold on 
plot(param.true.rho(:,ones(1,nem-1)),'linestyle','--','color',blue,'linewidth',1.5)
hold off
ylim([0,1])
xlim([0,nem])
xlabel('Iteration')   
ylabel('$\rho$')

subplot(4,2,7)
if isempty(strmatch('alpha',fixparam))
    set(gcf,'defaultaxescolororder',[blue;blue])
    plot(param.save.alpha(:,1:nem-1)','linewidth',1.5)
    hold on 
    plot(param.true.alpha(ones(1,nem-1),:),'--','linewidth',1.5)
    hold off
    xlim([0,nem])
    ylim([0,5])
    xlabel('Iteration')   
    ylabel('$\alpha$') 
elseif isempty(strmatch('beta',fixparam))
    plot(param.save.beta(:,1:nem-1)','linewidth',1)
    hold on 
    plot(param.true.beta(ones(1,nem-1),:),'--','linewidth',1)
    hold off
    xlim([0,nem])
    ylim([0,5])
    xlabel('Iteration')   
    ylabel('$\beta_c$') 
end
    

subplot(4,2,8)
plot(param.save.mu(1:nem-1),'color',blue,'linewidth',1.5)
hold on 
plot(param.true.mu(:,ones(1,nem-1)),'linestyle','--','color',blue,'linewidth',1.5)
hold off
ylim([-0.5,1])
xlim([0,nem])
xlabel('Iteration')   
ylabel('$\mu$')  
end