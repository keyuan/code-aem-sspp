% Reproduce figure 3 in the manuscript


%% D1

clear synthdata
load('./data/synthdb/data_synthdb_c_10_alpha_4.mat')

y1      = synthdata.y;
time1   = synthdata.time;
delta1  = synthdata.delta;
in1     = synthdata.in;
xtrue1  = synthdata.x;

%% D2

clear synthdata
load('./data/synthdb/data_synthdb_c_1_alpha_4.mat')

y2      = synthdata.y;
time2   = synthdata.time;
delta2  = synthdata.delta;
in2     = synthdata.in;
xtrue2  = synthdata.x;

%% D3

clear synthdata
load('./data/synthdb/data_synthdb_c_1_alpha_1.mat')

y3      = synthdata.y;
time3   = synthdata.time;
delta3  = synthdata.delta;
in3     = synthdata.in;
xtrue3  = synthdata.x;
%%

figure(1), clf

for jj = 1:3
    
    switch jj
        case 1, 
            y = y1; time = time1; in = in1; xtrue = xtrue1;
        case 2, 
            y = y2; time = time2; in = in2; xtrue = xtrue2;
        case 3, 
            y = y3; time = time3; in = in3; xtrue = xtrue3;
    end
    
    subplot(1,3,jj)
    idx = find(in==1);
    u = 0.5; base = 5;
    lower = base*ones(size(idx));
    upper = base +u*ones(size(idx));
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
    if jj == 1, ylabel('Data'); end
    hold off
end
set(gcf,'units','centimeters');
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2),17,10]);

matlabfrag('./fig/fig-3')


