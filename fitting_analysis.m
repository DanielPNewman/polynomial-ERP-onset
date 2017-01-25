%By Dan Newman ('DN'; dan.newman86@gmail.com)and Rafael Abe ('RA'; abe_rafael@hotmail.com) 
%30/01/2014. This script finds the optimal order level for the polynomial to describe each participant's
%CPP using a sum of deviations/order size trade-off analysis. 
%Algorithm adapted from "finding the best trade-off point on a curve" at
%http://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

clear all
close all
clc

% path ='S:\R-MNHS-SPP\Bellgrove-data\4. Dan Newman\Dots TOT study\'; %DN
path ='D:\Dots TOT study\'; %DN


subject_folder = {'TTN01','TTN02','TTN03','TTN04','TTN05','TTN06',...
    'TTN07','TTN08','TTN09','TTN10','TTN11','TTN12','TTN13','TTN14','TTN15','TTN16'...
    ,'TTN17','TTN18','TTN19','TTN20','TTN21','TTN22','TTN23'...
     'TTN24','TTN25','TTN26','TTN27','TTN28','TTN29','TTN30','TTN31','TTN32','TTN33','TTN34','TTN35',...
     'TTN36','TTN37','TTN38','TTN39','TTN40','TTN41','TTN42','TTN43', 'TTN44','TTN45','TTN46','TTN47',...
    'TTN48','TTN49','TTN50','TTN51','TTN52','TTN53','TTN54','TTN55','TTN56','TTN57','TTN58','TTN59'};

allsubj = {'TTN01_D','TTN02_D','TTN03_D','TTN04_D','TTN05_D','TTN06_D',...
    'TTN07_D','TTN08_D','TTN09_D','TTN10_D','TTN11_D','TTN12_D','TTN13_D','TTN14_D','TTN15_D','TTN16_D'...
    ,'TTN17_D','TTN18_D','TTN19_D','TTN20_D','TTN21_D','TTN22_D','TTN23_D'...
    ,'TTN24_D','TTN25_D','TTN26_D','TTN27_D','TTN28_D','TTN29_D','TTN30_D'...
    ,'TTN31_D','TTN32_D','TTN33_D','TTN34_D','TTN35_D',...
    'TTN36_D','TTN37_D','TTN38_D','TTN39_D','TTN40_D','TTN41_D','TTN42_D','TTN43_D',...
    'TTN44_D','TTN45_D','TTN46_D','TTN47_D','TTN48_D','TTN49_D','TTN50_D','TTN51_D','TTN52_D','TTN53_D','TTN54_D','TTN55_D','TTN56_D','TTN57_D','TTN58_D','TTN59_D'};
duds = [1,4,18,19,21,28,51]; % 21 missing .mat files. %28 missing eeg files. 10,12,16,18,19 have low number of valid trials - not enough for accurate time-on-task analysis
single_participants = []; %DN: put the number of the participants you want to include in here and it will only run the analysis on them

file_start = 1; % bear in mind duds.
allblocks = {[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],...
    [1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],...
    [1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],...
    [1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1]...
    [1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],...
    [1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1],[1:1]};


if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    allblocks([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    allblocks = allblocks(single_participants);
end

% 'Trigger 1: coherence 18, patch 1, motion dir 270, coh motion onset 1420
% 'Trigger 2: coherence 18, patch 2, motion dir 270, coh motion onset 1420
% 'Trigger 3: coherence 18, patch 1, motion dir 270, coh motion onset 1820
% 'Trigger 4: coherence 18, patch 2, motion dir 270, coh motion onset 1820
% 'Trigger 5: coherence 18, patch 1, motion dir 270, coh motion onset 2220
% 'Trigger 6: coherence 18, patch 2, motion dir 270, coh motion onset 2220

% side,iti
targcodes = zeros(2,3); %targcodes(side,iti)
targcodes(1,:) = [101 103 105]; % left patch,
targcodes(2,:) = [102 104 106]; % right patch


Fs=500;
numch=64;
rtlim=[0.2 2.2];

% ch = [25]; %Channel Pz actiCAP64
ch = [31]; %biosemi
ch_R_L = [64,60]; % right hemi channels for left target, vice versa.
ch_for_ipsicon(1,:) = [64;60];
ch_for_ipsicon(2,:) = [60;64];
ch_l = [60];
ch_r = [64];
% time scales:

% stim-locked erps
ts = -0.376*Fs:1.875*Fs;   % in sample points, the ERP epoch
t = ts*1000/Fs;

% response-locked erps:
% trs = [-1.125*Fs:Fs*.375];
% trs = [-.875*Fs:Fs*.500];
trs = [-.876*Fs:Fs*.125];
tr = trs*1000/Fs;

% T = [-100:50:1500]; % stim-locked stft timepoint in ms
% Tr = [-600:50:50];  % r-locked
fftlen = 253;   % number of samples. 253 is exactly 7 cycles of 85/6=14.167 Hz.
F = [0:25]*Fs/fftlen;

num_bins=6;
invalid=zeros(length(allsubj),num_bins);

avRT=[];
avERP = []; avERPr = []; %avLRP_L = []; avLRP_R = []; avLRPr_L = []; avLRPr_R = [];
for s=file_start:length(allsubj)
    erp=[]; erpr = []; trial=0; ERP_function_temp=[];trial2=0;
    
    load([path subject_folder{s} '\' allsubj{s} '_resample_erp_dan_rafa']); %DN
    load([path subject_folder{s} '\' allsubj{s} '_CPP_ROIs']); %Ger, you can just switch CPP ROI for a single electrode like Pz for now
    clear matfiles;
    matfile = [path subject_folder{s} '\' allsubj{s} '.mat'];    
    load(matfile);
    
    %     A=find(allTrig==107|allTrig==113|allTrig==119|allTrig==125|allTrig==131|allTrig==137|allTrig==143|allTrig==149|allTrig==155);
    %     allTrig(A)=101;
    %     B=find(allTrig==108|allTrig==114|allTrig==120|allTrig==126|allTrig==132|allTrig==138|allTrig==144|allTrig==150|allTrig==156);
    %     allTrig(B)=102;
    %     C=find(allTrig==109|allTrig==115|allTrig==121|allTrig==127|allTrig==133|allTrig==139|allTrig==145|allTrig==151|allTrig==157);
    %     allTrig(C)=103;
    %     D=find(allTrig==110|allTrig==116|allTrig==122|allTrig==128|allTrig==134|allTrig==140|allTrig==146|allTrig==152|allTrig==158);
    %     allTrig(D)=104;
    %     E=find(allTrig==111|allTrig==117|allTrig==123|allTrig==129|allTrig==135|allTrig==141|allTrig==147|allTrig==153|allTrig==159);
    %     allTrig(E)=105;
    %     F=find(allTrig==112|allTrig==118|allTrig==124|allTrig==130|allTrig==136|allTrig==142|allTrig==148|allTrig==154|allTrig==160);
    %     allTrig(F)=106;
    %     clear A B C D E F
    
    % DN: ignore the catch trials (155,156,157,158,159,160) for now
    A=find(allTrig==107|allTrig==113|allTrig==119|allTrig==125|allTrig==131|allTrig==137|allTrig==143|allTrig==149);  %Ger you'll just use normal allTrig vector here, this is just here because of the hack we used for catch trials (we made it so 90% of trials set at the participant's coherence level determined by the staircase procedure, and the other 10% set as catch trials)
    allTrig(A)=101;
    B=find(allTrig==108|allTrig==114|allTrig==120|allTrig==126|allTrig==132|allTrig==138|allTrig==144|allTrig==150);
    allTrig(B)=102;
    C=find(allTrig==109|allTrig==115|allTrig==121|allTrig==127|allTrig==133|allTrig==139|allTrig==145|allTrig==151);
    allTrig(C)=103;
    D=find(allTrig==110|allTrig==116|allTrig==122|allTrig==128|allTrig==134|allTrig==140|allTrig==146|allTrig==152);
    allTrig(D)=104;
    E=find(allTrig==111|allTrig==117|allTrig==123|allTrig==129|allTrig==135|allTrig==141|allTrig==147|allTrig==153);
    allTrig(E)=105;
    F=find(allTrig==112|allTrig==118|allTrig==124|allTrig==130|allTrig==136|allTrig==142|allTrig==148|allTrig==154);
    allTrig(F)=106;
    clear A B C D E F
    [B,A]=butter(4,8*2/Fs);
    disp('##################################');
    for trial1=1:size(erp,3)
         slope=0; r=0;
        if ~fixation_break(trial1) && allTrig(trial1)>0 %Ger, you probably won't have this fixation break vector so you can delete fixation_break here
            trial=trial+1;
            for order=1:30 %Ger the only thing we are unsure of is whether this should be set at 30 as it is or higher (say 40). This does effect the optimal order calculated, you can change it and run the script again and see what I mean
                if perf(trial)==1
                trial2=trial2+1;
                [P,S(trial2,order)]= polyfit(t,mean(erp(CPP_ROI,:,trial1),1),order); %Ger, you can just switch CPP ROI for electrode Pz for now
                norm(trial2,order)=S(trial2,order).normr;
                end
            end
            meannorm(s,:)=mean(norm,1);
        end
    end
     % curve = [8.4663 8.3457 5.4507 5.3275 4.8305 4.7895 4.6889 4.6833 4.6819 4.6542 4.6501 4.6287 4.6162 4.585 4.5535 4.5134 4.474 4.4089 4.3797 4.3494 4.3268 4.3218 4.3206 4.3206 4.3203 4.2975 4.2864 4.2821 4.2544 4.2288 4.2281 4.2265 4.2226 4.2206 4.2146 4.2144 4.2114 4.1923 4.19 4.1894 4.1785 4.178 4.1694 4.1694 4.1694 4.1556 4.1498 4.1498 4.1357 4.1222 4.1222 4.1217 4.1192 4.1178 4.1139 4.1135 4.1125 4.1035 4.1025 4.1023 4.0971 4.0969 4.0915 4.0915 4.0914 4.0836 4.0804 4.0803 4.0722 4.065 4.065 4.0649 4.0644 4.0637 4.0616 4.0616 4.061 4.0572 4.0563 4.056 4.0545 4.0545 4.0522 4.0519 4.0514 4.0484 4.0467 4.0463 4.0422 4.0392 4.0388 4.0385 4.0385 4.0383 4.038 4.0379 4.0375 4.0364 4.0353 4.0344];
    curve=(meannorm(s,:));
    %# get coordinates of all the points
    nPoints = length(curve);
    allCoord = [1:nPoints;curve]';              %'# SO formatting
    
    %# pull out first point
    firstPoint = allCoord(1,:);
    
    %# get vector between first and last point - this is the line
    lineVec = allCoord(end,:) - firstPoint;
    
    %# normalize the line vector
    lineVecN = lineVec / sqrt(sum(lineVec.^2));
    
    %# find the distance from each point to the line:
    %# vector between all points and first point
    vecFromFirst = bsxfun(@minus, allCoord, firstPoint);
    
    %# To calculate the distance to the line, we split vecFromFirst into two
    %# components, one that is parallel to the line and one that is perpendicular
    %# Then, we take the norm of the part that is perpendicular to the line and
    %# get the distance.
    %# We find the vector parallel to the line by projecting vecFromFirst onto
    %# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
    %# We project vecFromFirst by taking the scalar product of the vector with
    %# the unit vector that points in the direction of the line (this gives us
    %# the length of the projection of vecFromFirst onto the line). If we
    %# multiply the scalar product by the unit vector, we have vecFromFirstParallel
    scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
    vecFromFirstParallel = scalarProduct * lineVecN;
    vecToLine = vecFromFirst - vecFromFirstParallel;
    
    %# distance to line is the norm of vecToLine
    distToLine = sqrt(sum(vecToLine.^2,2));
    
    %# plot the distance to the line
    figure('Name','distance from curve to line'), plot(distToLine)
    
    %# now all you need is to find the maximum
    [maxDist,cpp_idxOfBestPoint] = max(distToLine);
    
    %# plot
    figure, plot(curve)
    hold on
    plot(allCoord(cpp_idxOfBestPoint,1), allCoord(cpp_idxOfBestPoint,2), 'or')
    
    subjID = [path subject_folder{s} '\' allsubj{s}];
    save([subjID '_CPP_order'],'cpp_idxOfBestPoint') %this saves the optimal order for each participant 
end
% curve = [8.4663 8.3457 5.4507 5.3275 4.8305 4.7895 4.6889 4.6833 4.6819 4.6542 4.6501 4.6287 4.6162 4.585 4.5535 4.5134 4.474 4.4089 4.3797 4.3494 4.3268 4.3218 4.3206 4.3206 4.3203 4.2975 4.2864 4.2821 4.2544 4.2288 4.2281 4.2265 4.2226 4.2206 4.2146 4.2144 4.2114 4.1923 4.19 4.1894 4.1785 4.178 4.1694 4.1694 4.1694 4.1556 4.1498 4.1498 4.1357 4.1222 4.1222 4.1217 4.1192 4.1178 4.1139 4.1135 4.1125 4.1035 4.1025 4.1023 4.0971 4.0969 4.0915 4.0915 4.0914 4.0836 4.0804 4.0803 4.0722 4.065 4.065 4.0649 4.0644 4.0637 4.0616 4.0616 4.061 4.0572 4.0563 4.056 4.0545 4.0545 4.0522 4.0519 4.0514 4.0484 4.0467 4.0463 4.0422 4.0392 4.0388 4.0385 4.0385 4.0383 4.038 4.0379 4.0375 4.0364 4.0353 4.0344];
curve=(mean(meannorm,1));
%# get coordinates of all the points
nPoints = length(curve);
allCoord = [1:nPoints;curve]';              %'# SO formatting

%# pull out first point
firstPoint = allCoord(1,:);

%# get vector between first and last point - this is the line
lineVec = allCoord(end,:) - firstPoint;

%# normalize the line vector
lineVecN = lineVec / sqrt(sum(lineVec.^2));

%# find the distance from each point to the line:
%# vector between all points and first point
vecFromFirst = bsxfun(@minus, allCoord, firstPoint);

%# To calculate the distance to the line, we split vecFromFirst into two 
%# components, one that is parallel to the line and one that is perpendicular 
%# Then, we take the norm of the part that is perpendicular to the line and 
%# get the distance.
%# We find the vector parallel to the line by projecting vecFromFirst onto 
%# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
%# We project vecFromFirst by taking the scalar product of the vector with 
%# the unit vector that points in the direction of the line (this gives us 
%# the length of the projection of vecFromFirst onto the line). If we 
%# multiply the scalar product by the unit vector, we have vecFromFirstParallel
scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
vecFromFirstParallel = scalarProduct * lineVecN;
vecToLine = vecFromFirst - vecFromFirstParallel;

%# distance to line is the norm of vecToLine
distToLine = sqrt(sum(vecToLine.^2,2));

%# plot the distance to the line
figure('Name','distance from curve to line'), plot(distToLine)

%# now all you need is to find the maximum
[maxDist,cpp_idxOfBestPoint] = max(distToLine);

%# plot
figure, plot(curve)
hold on
plot(allCoord(cpp_idxOfBestPoint,1), allCoord(cpp_idxOfBestPoint,2), 'or')


