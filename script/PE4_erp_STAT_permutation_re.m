%% load files
clc;clear
cd('G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult');

filename={'day1_nsubavg_re.mat','day2_nsubavg_re.mat','day3_nsubavg_re.mat','day4_nsubavg_re.mat'};
%%
for md=1:4
    load(filename{md});
    %reshape ????¡À???
    for isub=1:14
        for idev=1:2
    Diff{md,isub,idev}=DEV_avg{isub,idev};
        end
    end
end
%% statistics
% cfg_neighb        = [];
% cfg_neighb.method = 'distance';
% neighbours        = ft_prepare_neighbours(cfg_neighb,  STD2timelock4{1});
cfg = [];

layout = ft_prepare_layout(cfg,  STD_avg{1,1} );

cfg.method = 'triangulation';
% % cfg.neighbourdist = 0.7;
cfg.layout = layout;
cfg.feedback = 'yes'; %%with this you get the feedback plot
neighbours = ft_prepare_neighbours(cfg, STD_avg{1,1} );
%% P3

cfg=[];
cfg.channel     = {'F3', 'F4', 'Fz','C3', 'C4', 'Cz','P3', 'P4', 'Pz'};
%cfg.latency          = [0.25 0.35];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            =0.025/6;
cfg.numrandomization = 'all';
cfg.minnbchan        = 0;

Nsubj  = 14;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

for idev=1:2
    i=1;
    for iday1=1:3
        for iday2=iday1+1:4
            a=Diff(iday1,:,idev);b=Diff(iday2,:,idev);
            stat{i,idev}= ft_timelockstatistics(cfg,a{:}, b{:});
            i=i+1;
        end
    end
end
%%
cfg = [];
cfg.alpha  = 0.025/6;
cfg.parameter = 'stat';
cfg.zlim   = [-4 4];
cfg.layout = layout;
ft_clusterplot(cfg, stat);
for idev=1:2

    for i=1:6
 
           ft_clusterplot(cfg, stat{i,idev});

    end
end