% updated on 190717 to read multi-channel z stacks. The 488 bead images are
% Channel 3 in the 3-channel z stacks

% updated on 190314.

clear all
close all

% updated on 190221 to incorporate adaptive thresholding, and also use the
% latest version of fitMultipleFoci so that one can define the maximum
% number of foci identified in each FOV.

% update on 210125, I added this line in the middle: Xdrift = []
% if xdrift is not renewed, and break out of the for loop, the old drift 
% file will be save again for current FOV 


NFOV = 22; % number of fields of views
NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
InitialLocalMaxThresh = 100; % brightness threshold for bead identification
ImageSize = 1536; % number of pxls
AdjrsquareThreshold = 0.9;
NumHybs = 14;
MaxNumBeadsToFit = 50;
Hyb0IsBit1 = 0; % change this to 1 if hyb0 is bit1.
%%
if exist('ShiftList.mat')==2
    load('ShiftList.mat');
else
    ShiftList = zeros(NumHybs,2);
end
if ~exist('beadfigs', 'dir')
    mkdir('beadfigs');
end
if ~exist('TracingDriftParams', 'dir')
    mkdir('TracingDriftParams');
end

BadFOV = [];

for jj = 1:NFOV-1
    % skip the FOV that have problems
%     if ismember(jj,FOVtoExclude)
%         continue
%     end

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
    
    FileName = ['sequential/STORM1_00_' FOVid];
    [ImageStack, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,3); % updated on 190717
    Image1 = mean(ImageStack,3);
%     figure(2)
%     imagesc(Image1)
%     hold on
%     colormap gray
%     axis equal
    LocalMaxThresh = InitialLocalMaxThresh;
    [Xfit1, Yfit1, Zfit1, Xgof1, Ygof1, Zgof1, Intensity1, Xwidth, Ywidth, Zwidth] = fitMultipleFoci(ImageStack,LocalMaxThresh,MaxNumBeadsToFit);
    Ind = find([Xgof1.adjrsquare]>AdjrsquareThreshold & ...
            [Ygof1.adjrsquare]>AdjrsquareThreshold & ...
            [Zgof1.adjrsquare]>AdjrsquareThreshold);
    Xfit1 = Xfit1(Ind);
    Yfit1 = Yfit1(Ind);
    Zfit1 = Zfit1(Ind);
    NumLandMarks = length(Xfit1);
    while NumLandMarks<10 && LocalMaxThresh>15
        LocalMaxThresh = LocalMaxThresh-10;
        [Xfit1, Yfit1, Zfit1, Xgof1, Ygof1, Zgof1, Intensity1, Xwidth, Ywidth, Zwidth] = fitMultipleFoci(ImageStack,LocalMaxThresh,MaxNumBeadsToFit);
        Ind = find([Xgof1.adjrsquare]>AdjrsquareThreshold & ...
            [Ygof1.adjrsquare]>AdjrsquareThreshold & ...
            [Zgof1.adjrsquare]>AdjrsquareThreshold);
        Xfit1 = Xfit1(Ind);
        Yfit1 = Yfit1(Ind);
        Zfit1 = Zfit1(Ind);
        NumLandMarks = length(Xfit1);
    end
    display([num2str(NumLandMarks) ' beads identified.']);

    if NumLandMarks<2
        display('NumLandMarks<2')
        continue
    end
    
    % find the minimum 2D distance to each bead (from all other beads): MinDisArray
	MinDisMatrix = ones(NumLandMarks,NumLandMarks)*ImageSize;
	for i = 1:NumLandMarks
		for j = 1:NumLandMarks
            if i~=j
    			MinDisMatrix(i,j) = ((Xfit1(i)-Xfit1(j))^2+(Yfit1(i)-Yfit1(j))^2)^0.5;
            end
		end
    end
    
    % This contains the minimum value of each column.
	MinDisArray = min(MinDisMatrix);
    % update on 210125, if xdrift is not renewed, and break out of the
    % belowfor loop, the old drift file will be saved again for current FOV
    Xdrift = [];
    
    for ii = 1:NumHybs
        if Hyb0IsBit1 == 0
            if ii<10
                FileName = ['sequential/STORM1_0' num2str(ii) '_' FOVid];
            else 
                FileName = ['sequential/STORM1_' num2str(ii) '_' FOVid];
            end
        elseif Hyb0IsBit1 == 1
            if ii-1<10
                FileName = ['sequential/STORM1_0' num2str(ii-1) '_' FOVid];
            else
                FileName = ['sequential/STORM1_' num2str(ii-1) '_' FOVid];
            end
        end        
        [ImageStack, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,3); % updated on 190717
        Image2 = mean(ImageStack,3);
%         figure(3)
%         imagesc(Image2)
%         hold on
%         colormap gray
%         axis equal
        LocalMaxThresh = InitialLocalMaxThresh;
        [Xfit2, Yfit2, Zfit2, Xgof2, Ygof2, Zgof2, Intensity2, Xwidth, Ywidth, Zwidth] = fitMultipleFoci(ImageStack,LocalMaxThresh,MaxNumBeadsToFit);
        
        if isempty(Xfit2)
            display('!!!No beads found!!!')
            BadFOV = [BadFOV;FileName];
            break
        end 
        
        Ind = find([Xgof2.adjrsquare]>AdjrsquareThreshold & ...
            [Ygof2.adjrsquare]>AdjrsquareThreshold & ...
            [Zgof2.adjrsquare]>AdjrsquareThreshold);
        Xfit2 = Xfit2(Ind);
        Yfit2 = Yfit2(Ind);
        Zfit2 = Zfit2(Ind);
        NumLandMarks2 = length(Xfit2);
        while NumLandMarks2<10 && LocalMaxThresh>15
            LocalMaxThresh = LocalMaxThresh-10;
            [Xfit2, Yfit2, Zfit2, Xgof2, Ygof2, Zgof2, Intensity2, Xwidth, Ywidth, Zwidth] = fitMultipleFoci(ImageStack,LocalMaxThresh,MaxNumBeadsToFit);
            Ind = find([Xgof2.adjrsquare]>AdjrsquareThreshold & ...
                [Ygof2.adjrsquare]>AdjrsquareThreshold & ...
                [Zgof2.adjrsquare]>AdjrsquareThreshold);
            Xfit2 = Xfit2(Ind);
            Yfit2 = Yfit2(Ind);
            Zfit2 = Zfit2(Ind);
            NumLandMarks2 = length(Xfit2);
        end
        display([num2str(NumLandMarks2) ' beads identified.']);
        
%         figure(3)
%         plot(Xfit2,Yfit2,'x');
%         figure(2)
%         plot(Xfit2,Yfit2,'x');
        
        NumLandMarks2 = length(Xfit2);
        if NumLandMarks2<2
            display('NumLandMarks2<2')
            break
        end 
        
        % find the minimum 2D distance to each bead (from all other beads): MinDisArray
        MinDisMatrix2 = ones(NumLandMarks2,NumLandMarks2)*ImageSize;
        for i = 1:NumLandMarks2
            for j = 1:NumLandMarks2
                if i~=j
                    MinDisMatrix2(i,j) = ((Xfit2(i)-Xfit2(j))^2+(Yfit2(i)-Yfit2(j))^2)^0.5;
                end
            end
        end
        MinDisArray2 = min(MinDisMatrix2);
        
        %match molecules from Image1 and Image2 first in 2D and then in 3D
        Xfit1_match = [];
        Yfit1_match = [];
        Zfit1_match = [];
        Xfit2_match = [];
        Yfit2_match = [];
        Zfit2_match = [];
        N = 0;
        for i = 1:NumLandMarks
            Xfit1_new = Xfit1(i) + ShiftList(ii,1);
            Yfit1_new = Yfit1(i) + ShiftList(ii,2);
            [m, Ind] = min(((Xfit2-Xfit1_new).^2+(Yfit2-Yfit1_new).^2).^0.5);
            if ~isempty(m)
                if m<MinDisArray(i)/2 && m<MinDisArray2(Ind)/2 && m<60
                    N = N+1;
                    % The fitting for beads in hyb 0 is not changed, but
                    % the fitting in later hybs are those that have the
                    % shorted drifting.
                    Xfit1_match(N) = Xfit1(i);
                    Yfit1_match(N) = Yfit1(i);
                    Zfit1_match(N) = Zfit1(i);
                    Xfit2_match(N) = Xfit2(Ind);
                    Yfit2_match(N) = Yfit2(Ind);
                    Zfit2_match(N) = Zfit2(Ind);
                end
            end
        end
        
        if N>=2
            MeanDriftX = mean(Xfit2_match-Xfit1_match);
            MeanDriftY = mean(Yfit2_match-Yfit1_match);
            MeanDriftZ = mean(Zfit2_match-Zfit1_match);
            x_new = Xfit2_match-MeanDriftX;
            y_new = Yfit2_match-MeanDriftY;
            z_new = Zfit2_match-MeanDriftZ;
            MedianErrorX = median(Xfit1_match-x_new);
            MedianErrorY = median(Yfit1_match-y_new);
            MedianErrorZ = median(Zfit1_match-z_new);
            Ind = find(abs(Xfit1_match-x_new-MedianErrorX)<1 & abs(Yfit1_match-y_new-MedianErrorY)<1 ...
                & abs(Zfit1_match-z_new-MedianErrorZ)<1);
            if length(Ind)>=2
                Xfit1_match = Xfit1_match(Ind);
                Yfit1_match = Yfit1_match(Ind);
                Zfit1_match = Zfit1_match(Ind);
                Xfit2_match = Xfit2_match(Ind);
                Yfit2_match = Yfit2_match(Ind);
                Zfit2_match = Zfit2_match(Ind);
                MeanDriftX = mean(Xfit2_match-Xfit1_match);
                MeanDriftY = mean(Yfit2_match-Yfit1_match);
                MeanDriftZ = mean(Zfit2_match-Zfit1_match);
                x_new = Xfit2_match-MeanDriftX;
                y_new = Yfit2_match-MeanDriftY;
                z_new = Zfit2_match-MeanDriftZ;

                Ind = find(abs(Xfit1_match-x_new)<1 & abs(Yfit1_match-y_new)<1 & abs(Zfit1_match-z_new)<1);
                if length(Ind)>=2
                    Xfit1_match = Xfit1_match(Ind);
                    Yfit1_match = Yfit1_match(Ind);
                    Zfit1_match = Zfit1_match(Ind);
                    Xfit2_match = Xfit2_match(Ind);
                    Yfit2_match = Yfit2_match(Ind);
                    Zfit2_match = Zfit2_match(Ind);
                    MeanDriftX = mean(Xfit2_match-Xfit1_match);
                    MeanDriftY = mean(Yfit2_match-Yfit1_match);
                    MeanDriftZ = mean(Zfit2_match-Zfit1_match);
                    x_new = Xfit2_match-MeanDriftX;
                    y_new = Yfit2_match-MeanDriftY;
                    z_new = Zfit2_match-MeanDriftZ;
%                     figure(2)
%                     plot(Xfit1_match,Yfit1_match,'.');
%                     plot(x_new,y_new,'o');
                else
                    display('less than 2 beads with final error less than 1')
                    break
                end
            else 
                display('less than 2 beads left after restricting error of error to less than 1')
                break
            end  
        else
            display('less than 2 beads matched in 2D')
            break
        end

        Xdrift(ii) = MeanDriftX;
        Ydrift(ii) = MeanDriftY;
        Zdrift(ii) = MeanDriftZ;

    end
    % output
    if length(Xdrift) == NumHybs
        save(['TracingDriftParams\DriftParams' FOVid '.mat'], 'Xdrift', 'Ydrift', 'Zdrift');
%         figure(2)
%         savefig(['beadfigs/figs2_' FOVid '.fig']);
    end
end
save('BadFOV_noBeads.mat','BadFOV')
    



    