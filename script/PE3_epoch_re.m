clear;clc;

basedir = {'G:\passiveexp22222\data\4thanalysis\day1','G:\passiveexp22222\data\4thanalysis\day2',...
    'G:\passiveexp22222\data\4thanalysis\day3','G:\passiveexp22222\data\4thanalysis\day4'};
filt ='*_deIC3.set';
matfileNAM={'pre','post'};
deviant_types = {'2','3'};
standard_type = 11;
art_thresh = 200; % +/- uV
twin = [-0.1 0.6];
savedir='G:\passiveexp22222\data\4thanalysis\MVPA_dataNresult';
%%
eeglab;
close(gcf);
%%
for md=1:length(basedir)
    
    cd(basedir{md});
    savenam=strrep(basedir{md}, '\', '_');
    savenam=strrep(savenam, 'G:_passiveexp22222_data_4thanalysis_', '');
    outputdir = basedir{md};
    files = dir(filt);
    ICtotal_brain_good=[];
    %%
    for curfile =1:length(files)
        file = files(curfile).name;
        EEG = pop_loadset(file,pwd);
        [pth,nam,ext] = fileparts(file);
        
        nam = nam(1:19);
        fprintf('Working on %s\n',[nam ext]);
        %% reref
        % Remove ICA-related info (if present) and insert A2 at ch21.
        EEG.nbchan      = 21;
        EEG.data(21,:)  = zeros(1,EEG.pnts);
        EEG.icaact      = [];
        EEG.icawinv     = [];
        EEG.icasphere   = [];
        EEG.icaweights  = [];
        EEG.icachansind = [];
        EEG.chanlocs(21).label = 'A2';
        EEG = pop_chanedit(EEG, 'changefield',{21 'labels' 'A2'},'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc', ...
            'eval','chans = pop_chancenter( chans, [],[]);');
        EEG = pop_reref( EEG, [20 21] );
        %EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_ref4.set'],'filepath',outputdir);
        %EEG = pop_select( EEG, 'channel',{'C3','C4','Cz','F3','F4','P3','P4','Fz','Pz'});
        %%
        ind_boundary=[];
        for i=1:length(EEG.event)
            if ismember(str2num(EEG.event(i).type),[1,4,16,64])
                EEG.event(i).type='1';
                
            elseif ismember(str2num(EEG.event(i).type),[3,12,48,192])
                EEG.event(i).type='3';
                
            elseif ismember(str2num(EEG.event(i).type),[2,8,32,128])
                EEG.event(i).type='2';
            elseif ismember('boundary',EEG.event(i).type)
                ind_boundary=[ind_boundary;i];
            end
        end
        
        %% 去掉不必要的event
        
        for idev=1:length(deviant_types)
            %%
            uniqueToDEV =[];uniqueToSTD = [];indices1=[];indices2=[];indices3=[];indices4=[];
            [EEG1,indices0] = pop_epoch( EEG, { deviant_types{idev}  }, [-1         0.60]);
            
            
            
            [DEV,indices1] = pop_epoch( EEG1, { deviant_types{idev} }, [-0.1         0.60]);
            [STD,indices2] = pop_epoch( EEG1, { '1' }, [-0.1         0.60] );
            [DEV Indexes1]  = pop_eegthresh(DEV, 1,[3,4,9,10,11], -200, 200, DEV.xmin, DEV.xmax,1,1);
            [STD Indexes2]  = pop_eegthresh(STD, 1,[3,4,9,10,11], -200, 200, STD.xmin, STD.xmax,1,1);
            %EEG1 = pop_select( EEG1, 'channel',{'C3','C4','Cz','F3','F4','P3','P4','Fz','Pz'});
            
            %% remove baseline
            % [EEG1,indices] = pop_selectevent( EEG1, 'omitepoch',[uniqueToDEV;uniqueToSTD],'deleteevents','off','deleteepochs','on','invertepochs','off');
            DEV = pop_rmbase( DEV, [-100 0] ,[]);
            STD = pop_rmbase( STD, [-100 0] ,[]);
            ntrial.dev(curfile,idev)=DEV.trials;
            ntrial.std(curfile,idev)=STD.trials;
            
            DEV = eeglab2fieldtrip(DEV, 'preprocessing');
            STD = eeglab2fieldtrip(STD, 'preprocessing');
            
            %%
            cfg = [];
            cfg.keeptrials = 'yes';
            STD_keeptrials_nsub{curfile,idev}    = ft_timelockanalysis(cfg, STD);
            DEV_keeptrials_nsub{curfile,idev}    = ft_timelockanalysis(cfg, DEV);
            
            cfg = [];
            cfg.keeptrials = 'no';
            STD_avg{curfile,idev}   = ft_timelockanalysis(cfg, STD);
            DEV_avg{curfile,idev}   = ft_timelockanalysis(cfg, DEV);
            
            cfg           = [];
            cfg.operation = 'subtract';
            %             cfg.parameter = 'trial';
            %             clssi_Diff{curfile,idev} = ft_math(cfg, DEV_keeptrials_nsub{curfile,idev} , STD_keeptrials_nsub{curfile,idev} );
            cfg.parameter = 'avg';
            a=DEV_avg{curfile,idev} ;b=STD_avg{curfile,idev};
            Diff_avg{curfile,idev} = ft_math(cfg, a, b );
            %%

        end
        
    end
    save([savedir,'\',savenam,'_nsubavg_re'],'STD_avg','DEV_avg','Diff_avg','ntrial');
    save([savedir,'\',savenam,'_nsubDEVnSTDmvpa_re'],'STD_keeptrials_nsub','DEV_keeptrials_nsub');
    
end
% save([savedir,'\','peakdetection.mat'],'latency','amplitude');

% size(DEV_keeptrials_nsub{9,1}.trial)
% size(STD_keeptrials_nsub{9,1}.trial)