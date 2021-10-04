% updated on 190717 to read multi-channel z stacks.

% updated on 190314.

% updated on 190221: use adaptive thresholding to try to fit as high as MaxNumFociToFit
% foci per FOV. Later I can use the foci clustering information from all
% hybs to determine which foci are real and which foci are noise.

clear all
close all
NFOV =22; % number of fields of views
NumImage = 541; % Number of images in each dax
FramesToWait = 5; % frames to wait at each height for each channel
ImageSize = 1536; % number of pxls
AdjrsquareThreshold = 0.7;
NumHybs = 14; % number of secondary hybs
InitialLocalMaxThresh = 50; % brightness threshold for foci identification
MaxNumFociToFit = 50; % maxium number of foci to fit per FOV
FociAreaThreshold = 20; % this threshold eliminates large patches of fluorescence
TotalNumChannels = 3; % total number of channels in the multi-channel z stack

%%
if ~exist('TracingResults_SmallScale', 'dir')
    mkdir('TracingResults_SmallScale');
end
if ~exist('figs_FociFitting', 'dir')
    mkdir('figs_FociFitting');
end
load('DeltaZ.mat');
load('tform.mat');
for jj = 0:NFOV-1
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
    if exist(['TracingDriftParams/DriftParams' FOVid '.mat'])==2
        load(['TracingDriftParams/DriftParams' FOVid '.mat']);
    else
        continue
    end
    
    for i = 1:NumHybs
        if i<10
            FileName = ['sequential/STORM1_0' num2str(i) '_' FOVid];
        else
            FileName = ['sequential/STORM1_' num2str(i) '_' FOVid];
        end
        
        [ImageStack, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,1);     
                % Here warp the Cy5 image into the Cy3 channel.
        for j = 1:size(ImageStack,3)
            ImageStack(:,:,j) = imtransform(ImageStack(:,:,j), tform, 'XData', [1 ImageSize], 'Ydata', [1 ImageSize]);
        end
        % use adaptive thresholding to try to get MaxNumFociToFit foci
        LocalMaxThresh = InitialLocalMaxThresh;
        ImageMax = medfilt2(max(ImageStack, [],3));
        background = imopen(ImageMax, strel('disk', 4));
        ImageMax = ImageMax-background;
        ImageMax(find(ImageMax<0)) = 0;
        BW = imextendedmax(ImageMax,LocalMaxThresh);
        CC = bwconncomp(BW, 8);
        S = regionprops(CC, 'Area');
        Ind = find([S.Area]<FociAreaThreshold);
        S = S(Ind);
        while length(S)<MaxNumFociToFit && LocalMaxThresh>15
            LocalMaxThresh = LocalMaxThresh-10;
            BW = imextendedmax(ImageMax,LocalMaxThresh);
            CC = bwconncomp(BW, 8);
            S = regionprops(CC, 'Area');
            Ind = find([S.Area]<FociAreaThreshold);
            S = S(Ind);
        end
        
        [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity, Xwidth, Ywidth, Zwidth] = fitMultipleFoci(ImageStack,LocalMaxThresh,MaxNumFociToFit);
        if length(Xfit)>0
            Ind = find([Xgof.adjrsquare]>AdjrsquareThreshold & ...
                [Ygof.adjrsquare]>AdjrsquareThreshold & ...
                [Zgof.adjrsquare]>AdjrsquareThreshold);
            XfitList{i} = Xfit(Ind)- Xdrift(i);
            YfitList{i} = Yfit(Ind)- Ydrift(i);
            ZfitList{i} = Zfit(Ind)- Zdrift(i);
            ZfitList{i} = ZfitList{i}-DeltaZ; % cancel the color shift in Z
            XgofList{i} = Xgof(Ind);
            YgofList{i} = Ygof(Ind);
            ZgofList{i} = Zgof(Ind);
            IntensityList{i} = Intensity(Ind);
            XwidthList{i} = Xwidth(Ind);
            YwidthList{i} = Ywidth(Ind);
            ZwidthList{i} = Zwidth(Ind);
            figure(200)
            plot(XfitList{i},YfitList{i},'.');
        else
            XfitList{i} = [];
            YfitList{i} = [];
            ZfitList{i} = [];
            XgofList{i} = [];
            YgofList{i} = [];
            ZgofList{i} = [];
            IntensityList{i} = [];
            XwidthList{i} = [];
            YwidthList{i} = [];
            ZwidthList{i} = [];
        end
    end
    for i = NumHybs+1:NumHybs*2
        t = i-NumHybs;
        if t<10
            FileName = ['sequential/STORM1_0' num2str(t) '_' FOVid];
        else
            FileName = ['sequential/STORM1_' num2str(t) '_' FOVid];
        end
        
        [ImageStack, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,2); 
        % use adaptive thresholding to try to get MaxNumFociToFit foci
        LocalMaxThresh = InitialLocalMaxThresh;
        ImageMax = medfilt2(max(ImageStack, [],3));
        background = imopen(ImageMax, strel('disk', 4));
        ImageMax = ImageMax-background;
        ImageMax(find(ImageMax<0)) = 0;
        BW = imextendedmax(ImageMax,LocalMaxThresh);
        CC = bwconncomp(BW, 8);
        S = regionprops(CC, 'Area');
        Ind = find([S.Area]<FociAreaThreshold);
        S = S(Ind);
        while length(S)<MaxNumFociToFit && LocalMaxThresh>15
            LocalMaxThresh = LocalMaxThresh-10;
            BW = imextendedmax(ImageMax,LocalMaxThresh);
            CC = bwconncomp(BW, 8);
            S = regionprops(CC, 'Area');
            Ind = find([S.Area]<FociAreaThreshold);
            S = S(Ind);
        end
        
        [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity, Xwidth, Ywidth, Zwidth] = fitMultipleFoci(ImageStack,LocalMaxThresh,MaxNumFociToFit);
        if length(Xfit)>0
            Ind = find([Xgof.adjrsquare]>AdjrsquareThreshold & ...
                [Ygof.adjrsquare]>AdjrsquareThreshold & ...
                [Zgof.adjrsquare]>AdjrsquareThreshold);
            XfitList{i} = Xfit(Ind)- Xdrift(t);
            YfitList{i} = Yfit(Ind)- Ydrift(t);
            ZfitList{i} = Zfit(Ind)- Zdrift(t);
            XgofList{i} = Xgof(Ind);
            YgofList{i} = Ygof(Ind);
            ZgofList{i} = Zgof(Ind);
            IntensityList{i} = Intensity(Ind);
            XwidthList{i} = Xwidth(Ind);
            YwidthList{i} = Ywidth(Ind);
            ZwidthList{i} = Zwidth(Ind);
            figure(200)
            plot(XfitList{i},YfitList{i},'.');
        else
            XfitList{i} = [];
            YfitList{i} = [];
            ZfitList{i} = [];
            XgofList{i} = [];
            YgofList{i} = [];
            ZgofList{i} = [];
            IntensityList{i} = [];
            XwidthList{i} = [];
            YwidthList{i} = [];
            ZwidthList{i} = [];
        end
    end
    display([newline '*** Successfully fit ' num2str(length(XfitList{i})) ' foci in FOV_' FOVid ' ***' newline])
    figure(200)
    hold off
    savefig(['figs_FociFitting/FociFitting_' FOVid '.fig']);
    save(['TracingResults_SmallScale/result' FOVid '.mat'], 'XfitList', 'YfitList', 'ZfitList', 'XgofList', 'YgofList', 'ZgofList','IntensityList','XwidthList','YwidthList','ZwidthList');
end
