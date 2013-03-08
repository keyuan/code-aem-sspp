function plot_data(data,param)

time  = data.time;
xtrue = data.x;
y     = data.y;
in    = data.in;

set(0,'defaulttextinterpreter','latex');

red = [0.9, 0, 0];
blue = [0, 0, 0.7];
green = [0, 0.7, 0];

figure(1), clf

subplot(3,4,1:8)

idx = find(in==1);
u = 0.5;
lower = 5*ones(size(idx));
upper = 5+u*ones(size(idx));
plot([time(idx);time(idx)],[lower;upper],'color',blue)
hold on

plot(time,xtrue,'-.','color',red)
hold on

for ii = 1:size(y,1)
    idx = find(y(ii,:)==1);
    u = 0.2;
    lower = -1-(ii)*u*ones(size(idx));
    upper = -1-(ii+1)*u*ones(size(idx));
    plot([time(idx);time(idx)],[lower;upper],'color','k')
    hold on
end
ylim([-4,6])
xlabel('Time (s)')
ylabel('$\lambda_k\Delta$')

hold off

%= free energy convergence performance
subplot(3,4,9)
plot(lbsave(2:nem),'color',blue,'linewidth',1.5)

xlabel('Iteration')   
ylabel('$\mathcal{Q}(\theta,q)$')

%= parameter convergence performance
subplot(3,4,10)
plot(param.save.rho(1:nem-1),'color',blue,'linewidth',1.5)
hold on 
plot(param.true.rho(:,ones(1,nem-1)),'linestyle','--','color',blue,'linewidth',1.5)
hold off
ylim([0,1])
xlabel('Iteration')   
ylabel('$\rho$')

subplot(3,4,11)
set(gcf,'defaultaxescolororder',[blue;blue])
plot(param.save.alpha(:,1:nem-1)','linewidth',1.5)
hold on 
plot(param.true.alpha(ones(1,nem-1),:),'--','linewidth',1.5)
hold off
ylim([0,5])
xlabel('Iteration')   
ylabel('$\alpha$') 

subplot(3,4,12)
plot(param.save.mu(1:nem-1),'color',blue,'linewidth',1.5)
hold on 
plot(param.true.mu(:,ones(1,nem-1)),'linestyle','--','color',blue,'linewidth',1.5)
hold off
ylim([-1,1])
xlabel('Iteration')   
ylabel('$\mu$')  
end