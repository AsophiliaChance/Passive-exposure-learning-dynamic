% resample,high pass low pass and remove line noise,change channel
% name,channel location
clear;clc;
%%
resamp_srate = 250;
datadir = {'E:\passiveexp22222\data\4thanalysis\day1\TEST','E:\passiveexp22222\data\4thanalysis\day2',...
    'E:\passiveexp22222\data\4thanalysis\day3','E:\passiveexp22222\data\4thanalysis\day4'};
filt ='*.set';
deviant_types = {'2','3'};
%%
eeglab;
close(gcf);

%%
for md=1%:length(datadir)
    cd(datadir{md});
    outputdir =pwd;
    files = dir(filt);
    %%
    for curfile = 1:length(files)
        cd(datadir{md});
        file = files(curfile).name;
        EEG = pop_loadset(file,pwd);
        [pth,nam,ext] = fileparts(file);
        
        filenam = nam;
        fprintf('Working on %s\n',[nam ext]);
        
        
        EEG = pop_resample( EEG, resamp_srate);
        
        ft_data=eeglab2fieldtrip(EEG ,'preprocessing', 'none');
        
            cfg.lpfilter        = 'yes';
            cfg.lpfreq          =  30; 
            cfg.lpfiltdir  = 'twopass';
            cfg.lpfiltord     = 4;
        %     cfg.hpfilter        = 'yes';           
        %     cfg.hpfreq          =  0.1;           
        %     cfg.hpfiltdir  = 'twopass';
        %     cfg.hpfiltord     = 4;
        cfg.detrend     = 'yes';
        cfg.bsfilter='yes';
        cfg.bsfreq=[49 51;99 101];
        tempdata = ft_preprocessing(cfg,ft_data);
        
        %eeg_data= fieldtrip2eeglab(tempdata.hdr,cat(3,tempdata.trial{:}));
        EEG.data=tempdata.trial{1,1};
        %EEG = pop_select( EEG, 'nochannel',21);
        %EEG = pop_select( EEG, 'nochannel',20);
        for ichan=1:EEG.nbchan
            EEG.chanlocs(ichan).labels=EEG.chanlocs(ichan).labels(3:end);
        end
        
        EEG=pop_chanedit(EEG, 'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
        EEG = pop_saveset( EEG, 'filename',[filenam,'_nodelete1.set'],'filepath',outputdir);
        
        %% anti blank
        temp = struct2cell(EEG.event.').'; type = temp(:, 7); clear temp;%find(ismember(type,'boundary'));
        ind=find(ismember(type,'0')==0);ind=[1;ind];
        latency = [EEG.event.latency].';
        
        latenci=latency(ind);
        latenci_a=latenci;latenci_a(end)=[];
        latenci_b=latenci;latenci_b(1)=[];
        tmp=latenci_b-latenci_a;
        %temp=sort(tmp,'descend');
        delect_ind=find(tmp>195);
        
        pre_latenci=[latenci(delect_ind);latenci(end)+200];
        post_latenci=[latenci(delect_ind+1);length(EEG.times)];
        pre_latenci(:,2)=post_latenci;
        
        EEG = eeg_eegrej( EEG, pre_latenci);%934351
        
        %% Step5: Detect artefact
        %%automated FFT based continuous data rejection
        [~,rejections,~] = pop_rejcont(EEG,'onlyreturnselection','on','taper','hamming');
        rejdataPRO{md,curfile}=rejections;
        %                     toplot = [rejections repmat([0 1 0],size(rejections,1),1) repmat(1:21,size(rejections,1),1)];
        %                     eegplot(EEG.data,'srate',EEG.srate,'winlength',5,'winrej',toplot,'eloc_file',EEG.chanlocs);
        xmax(md,curfile,1)=EEG.xmax;
        EEG = pop_select(EEG,'nopoint',rejections);
        xmax(md,curfile,2)=EEG.xmax;
        xmax(md,curfile,3)=xmax(md,curfile,1)-xmax(md,curfile,2);
        
        %% 删除无用电极
        EEG = pop_select( EEG, 'nochannel',{'EOG'});
        %% 改mark
        
        for i=1:length(EEG.event)
            if ismember(str2num(EEG.event(i).type),[1,4,16,64])
                EEG.event(i).type='1';
                
            elseif ismember(str2num(EEG.event(i).type),[3,12,48,192])
                EEG.event(i).type='3';
                
            elseif ismember(str2num(EEG.event(i).type),[2,8,32,128])
                EEG.event(i).type='2';
            end
        end
        
        %% Step 10: Run ICA
        EEG_for_ICA = pop_eegfiltnew(EEG, 'locutoff',1);
        %%
        % EEG = pop_epoch( EEG, {  '2'  '3'  }, [-1         0.6], 'newname', ' resampled pruned with ICA epochs', 'epochinfo', 'yes');
        % EEG = pop_rmbase( EEG, [-100 0] ,[]);
        
        EEG_forICA = pop_resample(EEG_for_ICA, 100);
        
        EEG_forICA = pop_runica(EEG_forICA, 'icatype', 'runica', 'extended',1,'interrupt','off');
        EEG.icaweights = EEG_forICA.icaweights;
        EEG.icasphere  = EEG_forICA.icasphere;
        EEG = pop_saveset( EEG, 'filename',[filenam,'_ICA2.set'],'filepath',outputdir);
        
        
        
        
    end
end

