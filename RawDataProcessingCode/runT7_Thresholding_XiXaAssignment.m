clear all
close all
TotalNumTADs = 28;
ImageSize = 1536;
NFOV = 22; % number of fields of views
UmPerPxl = 0.108;
StepSize = 0.2; %um
NumImage = 541; % Number of images in each dax
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack

load('DeltaZ.mat');
load('tform.mat');
if ~exist('ChrChecking', 'dir')
    mkdir('ChrChecking');
end

%% analyze mean spatial distance matrix
n = 0;
for jj = 0:NFOV-1
% for jj = 2
    if NFOV<=10
        FOVid = ['0' num2str(jj)];
    elseif NFOV>10 && NFOV<=100
        if jj<10
            FOVid = ['0' num2str(jj)];
        else
            FOVid = [num2str(jj)];
        end
    elseif NFOV>100
        if jj<10
            FOVid = ['00' num2str(jj)];
        elseif jj<100
            FOVid = ['0' num2str(jj)];
        else
            FOVid = [num2str(jj)];
        end
    end 
    if exist(['Traces_SmallScale/TraceArrayRefined' FOVid '.mat'])==2
        load(['Traces_SmallScale/TraceArrayRefined' FOVid '.mat']);
    else
        continue
    end
    FileName = ['sequential/STORM1_00_' FOVid];
    [ImageStack_cy5, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,1);
    for j = 1:size(ImageStack_cy5,3)
        ImageStack_cy5(:,:,j) = imtransform(ImageStack_cy5(:,:,j), tform, 'XData', [1 ImageSize], 'Ydata', [1 ImageSize]);
    end
    FileName = ['sequential/STORM1_00_' FOVid];
    [ImageStack_cy3, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,2);
    
    
%     clear Chr
%     n = 0;
    ImageMax_cy3 = max(ImageStack_cy3,[],3);
%     ImageMax_cy5 = max(ImageStack_cy5,[],3);
    ImageMax_cy5 = mean(ImageStack_cy5,3);
    figure(1)
    f(1) = subplot(2,2,1);
    imagesc(ImageMax_cy5)
    colormap gray
    axis square
    f(2) = subplot(2,2,2);
    imagesc(ImageMax_cy3)
    colormap gray
    axis square
    f(3) = subplot(2,2,3);
    imagesc(ImageMax_cy5)
    colormap gray
    axis square
    hold on
    for i = 1:length(TraceArray)
        n = n+1;
        Cy5_value_list = [];
        Cy3_value_list = [];
        Chr(n).x = zeros(TotalNumTADs,1);
        Chr(n).y = zeros(TotalNumTADs,1);
        Chr(n).z = zeros(TotalNumTADs,1);
        Chr(n).r = zeros(TotalNumTADs,1);
        Chr(n).x(TraceArray{i}(:,end)) = TraceArray{i}(:,1);
        Chr(n).y(TraceArray{i}(:,end)) = TraceArray{i}(:,2);
        Chr(n).z(TraceArray{i}(:,end)) = TraceArray{i}(:,3);
        Chr(n).r(TraceArray{i}(:,end)) = 1;
        for k = 1:length(TraceArray{i})
            CurrentTrace = TraceArray{i};
            TADIdx = CurrentTrace(k,11);
            x = round(CurrentTrace(k,1)/UmPerPxl);
            y = round(CurrentTrace(k,2)/UmPerPxl);
            z = round(CurrentTrace(k,3)/StepSize);
            if z > 36 % some foci has fitting z = 37
                z = 36;
            end
            if x > 1536 
                x = 1536;
            end
            plot(x,y,'x');
            Cy5_value_list = [Cy5_value_list ImageStack_cy5(y,x,z)];
            Cy3_value_list = [Cy3_value_list ImageStack_cy3(y,x,z)];
        end
        Cy5_mean = mean(Cy5_value_list);
        Chr(n).cy5List = Cy5_value_list;
        Chr(n).cy5Mean = Cy5_mean;
        Cy3_mean = mean(Cy3_value_list);
        Chr(n).cy3List = Cy3_value_list;
        Chr(n).cy3Mean = Cy3_mean;
        text(x,y,[num2str(i) ',' num2str(round(Cy3_mean)) ',' num2str(round(Cy5_mean))],'Color','y');
    end
    set(f(1),'position',[0 .5 .5 .5])
    set(f(2),'position',[0 0 .5 .5])
    set(f(3),'position',[0.2 0 1 1])
    figure(1)
    hold off
    savefig(['ChrChecking/Probes_traces_' FOVid '.fig']);
end

%%
Xist_list = [];
Probe_list = [];
for k = 1:length(Chr)
    Xist_list = [Xist_list Chr(k).cy5Mean];
    Probe_list = [Probe_list Chr(k).cy3Mean];
end
%%
figure(11)
subplot(1,2,1)
[N,X] = hist(Probe_list,100);
str = '#46A928';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
Bh = bar(X,N,'facecolor',color);
hold on
title('Mean probe signal, cy3');
subplot(1,2,2)
[N,X] = hist(Xist_list,100);
str = '#EE5264';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
Bh = bar(X,N,'facecolor',color);
title('Mean Xist signal, cy5');
hold off
savefig(['hist plot of signal value']);
%%
XiIdx = [];
XaIdx = [];
XXIdx = [];
% set threshold according to the hist plot above
for k = 1:length(Chr)
    if Chr(k).cy3Mean > 1000 && Chr(k).cy5Mean < 500
        XaIdx = [XaIdx k];
    elseif Chr(k).cy3Mean > 1000 && Chr(k).cy5Mean > 800
        XiIdx = [XiIdx k];
    elseif Chr(k).cy3Mean > 1000
        XXIdx = [XXIdx k];
    end
end
AllXiChr = Chr(XiIdx);
AllXaChr = Chr(XaIdx);
AllXXChr = Chr(XXIdx);
save('AllXiChr','AllXiChr')
save('AllXaChr','AllXaChr')
save('AllXXChr','AllXXChr')
% save('Chr_BeforeXiAssignment','Chr')
display(['***   ' num2str(length(XiIdx)) ' Xi traces finally saved  ***'])
display(['***   ' num2str(length(XaIdx)) ' Xa traces finally saved  ***'])
display(['***   ' num2str(length(XXIdx)) ' XX traces finally saved  ***'])

