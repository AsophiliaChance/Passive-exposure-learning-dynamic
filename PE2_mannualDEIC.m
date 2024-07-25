%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DAY1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=1
EEG = pop_loadset('filename','2014_10_06T171430_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [7], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_10_06T171430_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  6], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_10_06T171430_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  3  8], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_10_20T171239_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  4], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_10_20T171239_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_10_20T171239_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  4  20], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');

%%
EEG = pop_loadset('filename','2014_11_24T170647_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  3  14], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_11_24T170647_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_11_24T170647_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_11_24T170647_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_12_08T170049_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_12_08T170049_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%% 13
EEG = pop_loadset('filename','2014_12_08T170049_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
%%
EEG = pop_loadset('filename','2014_12_08T170049_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
EEG = pop_subcomp( EEG, [2  5], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day1\\');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DAY2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2

EEG = pop_loadset('filename','2014_10_07T165604_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  7], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_10_07T165604_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  8], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_10_07T165604_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_10_21T171423_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  5], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%% 5
EEG = pop_loadset('filename','2014_10_21T171423_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  4], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_10_21T171423_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  16], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_11_25T165955_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [3  11], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_11_25T165955_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_11_25T165955_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [2  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');

%%
EEG = pop_loadset('filename','2014_11_25T165955_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_12_09T165207_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_12_09T165207_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_12_09T165207_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  4], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
%%
EEG = pop_loadset('filename','2014_12_09T165207_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
EEG = pop_subcomp( EEG, [1  2  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day2\\');
end
%%                               DAY3
for i=3
EEG = pop_loadset('filename','2014_10_08T164218_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  6], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 2
EEG = pop_loadset('filename','2014_10_08T164218_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  8], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 3
EEG = pop_loadset('filename','2014_10_08T164218_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  9], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 4
EEG = pop_loadset('filename','2014_10_22T165956_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  5], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 5
EEG = pop_loadset('filename','2014_10_22T165956_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 6
EEG = pop_loadset('filename','2014_10_22T165956_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  7], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 7
EEG = pop_loadset('filename','2014_11_26T170430_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  4], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 8
EEG = pop_loadset('filename','2014_11_26T170430_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 9
EEG = pop_loadset('filename','2014_11_26T170430_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 10
EEG = pop_loadset('filename','2014_11_26T170430_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 11
EEG = pop_loadset('filename','2014_12_10T164713_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 12
EEG = pop_loadset('filename','2014_12_10T164713_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 13
EEG = pop_loadset('filename','2014_12_10T164713_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  4], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
%% 14
EEG = pop_loadset('filename','2014_12_10T164713_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day3\\');
end
%%                               DAY4
for i=4
EEG = pop_loadset('filename','2014_10_09T165510_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  9], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 2
EEG = pop_loadset('filename','2014_10_09T165510_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  6], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 3
EEG = pop_loadset('filename','2014_10_09T165510_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 4
EEG = pop_loadset('filename','2014_10_23T164631_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  7], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 5
EEG = pop_loadset('filename','2014_10_23T164631_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 6
EEG = pop_loadset('filename','2014_10_23T164631_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  7], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%%  7
EEG = pop_loadset('filename','2014_11_27T165811_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  7], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 8
EEG = pop_loadset('filename','2014_11_27T165811_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  3], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 9
EEG = pop_loadset('filename','2014_11_27T165811_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 10
EEG = pop_loadset('filename','2014_11_27T165811_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  4  5], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 11
EEG = pop_loadset('filename','2014_12_11T165709_1_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [2  5], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 12
EEG = pop_loadset('filename','2014_12_11T165709_2_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  4], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 13
EEG = pop_loadset('filename','2014_12_11T165709_3_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  5], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
%% 14
EEG = pop_loadset('filename','2014_12_11T165709_4_ICA2.set','filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
EEG = pop_subcomp( EEG, [1  2], 0);
EEG = pop_saveset( EEG, 'filename',[EEG.filename(1:19),'_deIC3.set'],'filepath','E:\\passiveexp22222\\data\\4thanalysis\\day4\\');
end