% updated on 20201104, this program select Xi and Xa chr by manully cropingaccording 
% to Xist signal, primary probe signal and traces position

clear all
close all

NFOV = 22; % number of fields of views
NumImage = 541; % Number of images in each dax
FramesToWait = 5; % frames to wait at each height for each channel
ImageSize = 1536; % number of pxls
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
TotalNumTADs = 28;
UmPerPxl = 0.108;
StepSize = 0.2; %um
%%
% create the needed directoies
if ~exist('roiList_Xa', 'dir')
    mkdir('roiList_Xa');
end
if ~exist('roiList_Xi', 'dir')
    mkdir('roiList_Xi');
end

load('DeltaZ.mat');
load('tform.mat');
XiNum = 0;
XaNum = 0;
XXNum = 0;
for jj =1:NFOV-1
    if jj-1<10
        SampleNum = ['0' num2str(jj)];
    else
        SampleNum = num2str(jj);
    end
    if exist(['Traces_SmallScale/TraceArrayRefined' SampleNum '.mat']);
        load(['TracingDriftParams\DriftParams' SampleNum '.mat']);
        load(['Traces_SmallScale/TraceArrayRefined' SampleNum '.mat']);
    else
        continue
    end
    FileName = ['sequential/STORM1_00_' SampleNum];
    [ImageStack_cy5, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,1);
    for j = 1:size(ImageStack_cy5,3)
        ImageStack_cy5(:,:,j) = imtransform(ImageStack_cy5(:,:,j), tform, 'XData', [1 ImageSize], 'Ydata', [1 ImageSize]);
    end
    [ImageStack_cy3, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,2);
    
    ImageMax_cy5 = max(ImageStack_cy5,[],3);
    ImageMax_cy3 = max(ImageStack_cy3,[],3);
%     ImageMax_cy5 = mean(ImageStack_cy5,3);
%     ImageMax_cy3 = mean(ImageStack_cy3,3);
    %%
%     figure(99)
%     imagesc(ImageMax_cy3)
%     colormap gray
%     axis square
%     hold on
    % plot Xist signal on the left
    figure(100)
    f(1) = subplot(2,2,1);
    imagesc(ImageMax_cy5)
    colormap gray
    axis square
    hold on
    % add grid to the figure, easy to compare the foci position
    for row = 150:150:1536
        line([1, 1536], [row, row], 'Color', 'r', 'LineStyle',':','LineWidth',0.1);
    end
    for col  = 150:150:1536
        line([col, col], [1, 1536], 'Color', 'r', 'LineStyle',':','LineWidth',0.1);
    end
    hold off
    f(2) = subplot(2,2,3);
    imagesc(ImageMax_cy3)
    colormap gray
    axis square
    hold on
    for row = 150:150:1536
        line([1, 1536], [row, row], 'Color', 'r', 'LineStyle',':','LineWidth',0.1);
    end
    for col  = 150:150:1536
        line([col, col], [1, 1536], 'Color', 'r', 'LineStyle',':','LineWidth',0.1);
    end
    hold off
    
    %%
    % plot probe signal and traces on the right
    n = 0;
    figure(100)
    f(3) = subplot(2,2,2);
    imagesc(ImageMax_cy3)
    colormap gray
    axis square
    hold on
    % add grid to the figure, easy to compare the foci position
    for row = 150:150:1536
        line([1, 1536], [row, row], 'Color', 'r', 'LineStyle',':','LineWidth',0.1);
    end
    for col  = 150:150:1536
        line([col, col], [1, 1536], 'Color', 'r', 'LineStyle',':','LineWidth',0.1);
    end
    for i = 1:length(TraceArray)
        n = n+1;
        X_list = [];
        Y_list = [];
        %             Cy5_value_list = [];
        %             Cy3_value_list = [];
        Chr(n).x = zeros(TotalNumTADs,1);
        Chr(n).y = zeros(TotalNumTADs,1);
        Chr(n).z = zeros(TotalNumTADs,1);
        Chr(n).r = zeros(TotalNumTADs,1);
        Chr(n).x(TraceArray{i}(:,end)) = TraceArray{i}(:,1);
        Chr(n).y(TraceArray{i}(:,end)) = TraceArray{i}(:,2);
        Chr(n).z(TraceArray{i}(:,end)) = TraceArray{i}(:,3);
        Chr(n).r(TraceArray{i}(:,end)) = 1;
        for k = 1:size(TraceArray{i},1)
            CurrentTrace = TraceArray{i};
            TADIdx = CurrentTrace(k,11);
            x = round(CurrentTrace(k,1)/UmPerPxl);
            y = round(CurrentTrace(k,2)/UmPerPxl);
            X_list = [X_list x];
            Y_list = [Y_list y];
            %                 z = round(CurrentTrace(k,3)/StepSize);
            %                 if z > 36
            %                     continue
            %                 end
            plot(x,y,'.');
            %                 Cy5_value_list = [Cy5_value_list ImageStack_cy5(y,x,z)];
            %                 Cy3_value_list = [Cy3_value_list ImageStack_cy3(y,x,z)];
        end
        Chr(n).X_list = X_list;
        Chr(n).Y_list = Y_list;
        %             Cy5_mean = mean(Cy5_value_list);
        %             Chr(n).cy5List = Cy5_value_list;
        %             Chr(n).cy5Mean = Cy5_mean;
        %             Cy3_mean = mean(Cy3_value_list);
        %             Chr(n).cy3List = Cy3_value_list;
        %             Chr(n).cy3Mean = Cy3_mean;
        %             text(x,y,[num2str(i) ',' num2str(round(Cy3_mean)) ',' num2str(round(Cy5_mean))],'Color','y');
    end
    set(f(1),'position',[0 .5 .5 .5])
    set(f(2),'position',[0 0 .5 .5])
    set(f(3),'position',[0.2 0 1 1])

    % select Xi and Xa traces
    clear XiroiList
    clear XaroiList
    NewROI = questdlg('Do you want to select new ROI for Xi?'); %ask question
    if(strcmp(NewROI, 'Yes'))
        selectROI
        XiroiList = roiList;
        save(['roiList_Xi/roiList' SampleNum '.mat'],'roiList');
    else
        load(['roiList_Xi/roiList' SampleNum '.mat'])
        XiroiList = roiList;
    end
    NewROI = questdlg('Do you want to select new ROI for Xa?'); %ask question
    if(strcmp(NewROI, 'Yes'))
        selectROI
        XaroiList = roiList;
        save(['roiList_Xa/roiList' SampleNum '.mat'],'roiList');
    else
        load(['roiList_Xa/roiList' SampleNum '.mat'])
        XaroiList = roiList;
    end
    %% find the traces in Xi and Xa region
    % the the pixel list of croped Xi and Xa region
    CropPixelList_Xi = {};
    for i = 1:length(XiroiList)
        roiNOW = XiroiList(i);
        x1 = round(roiNOW.rect(1));
        x2 = round(roiNOW.rect(1)+roiNOW.rect(3));
        y1 = round(roiNOW.rect(2));
        y2 = round(roiNOW.rect(2)+roiNOW.rect(4));
        CurrPixelList = [];
        for x = x1:x2
            for y = y1:y2
                CurrPixel = 1536*(x-1)+y;
                CurrPixelList = [CurrPixelList CurrPixel];
            end
        end
        CropPixelList_Xi{i} = CurrPixelList;
    end
    CropPixelList_Xa = {};
    for i = 1:length(XaroiList)
        roiNOW = XaroiList(i);
        x1 = round(roiNOW.rect(1));
        x2 = round(roiNOW.rect(1)+roiNOW.rect(3));
        y1 = round(roiNOW.rect(2));
        y2 = round(roiNOW.rect(2)+roiNOW.rect(4));
        CurrPixelList = [];
        for x = x1:x2
            for y = y1:y2
                CurrPixel = 1536*(x-1)+y;
                CurrPixelList = [CurrPixelList CurrPixel];
            end
        end
        CropPixelList_Xa{i} = CurrPixelList;
    end
    %%
    % compare pixel list of croped region and traces, assignment
    for i = 1:length(Chr)
        XiOverlap = [];
        XaOverlap = [];
        
        CurrTracePixelList = (Chr(i).X_list-1)*1536 + Chr(i).Y_list;
        for j = 1:length(CropPixelList_Xi)
            XiOverlap = [XiOverlap length(intersect(CurrTracePixelList,CropPixelList_Xi{1,j}))];
        end
        for j = 1:length(CropPixelList_Xa)
            XaOverlap = [XaOverlap length(intersect(CurrTracePixelList,CropPixelList_Xa{1,j}))];
        end
        if max(XiOverlap)> max(XaOverlap)
            XiNum = XiNum + 1;
            AllXiChr(XiNum) = Chr(i);
        elseif max(XiOverlap)< max(XaOverlap)
            XaNum = XaNum + 1;
            AllXaChr(XaNum) = Chr(i);
        else
            XXNum = XXNum + 1;
            AllXXChr(XXNum) = Chr(i);
        end
    end
end

save('AllXiChr','AllXiChr')
save('AllXaChr','AllXaChr')
save('AllXXChr','AllXXChr')
display(['***   ' num2str(XiNum) ' Xi traces finally saved  ***'])
display(['***   ' num2str(XaNum) ' Xa traces finally saved  ***'])


































