% calibrate in 3D the two channels. Select 20 beads.
% generate the warping matrix from Cy5 channel (right) to the Cy3 channel (left).

clear all
close all

NumImage = 541; % Number of images in each dax
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
InitialLocalMaxThresh = 100; % brightness threshold for bead identification
ImageSize = 1536; % number of pxls
AdjrsquareThreshold = 0.9;
InitialLocalMaxThresh = 100; % Initial brightness threshold for bead identification

%%
lmol_match =[];
rmol_match =[];

for jj = 1:5
    %%
    FOVid = num2str(jj);
    FileName = ['calib/movie_000' FOVid];
    % 561 left; 647 right
    [ImageStack_left, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,2);
    [ImageStack_right, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName,NumImage,FramesToWait,TotalNumChannels,1);
    
    ImageMeanLeft = mean(ImageStack_left,3);
    ImageMeanRight = mean(ImageStack_right,3);
    
    %The bead intensity is normalized to their own maximum. Thus it doesn't
    % matter if the left and right channels have different signal intensity.
    RGB(:,:,1) = ImageMeanLeft/max(max(ImageMeanLeft));
    RGB(:,:,2) = ImageMeanRight/max(max(ImageMeanRight));
    RGB(:,:,3) = zeros(ImageSize,ImageSize);
    RGB(find(RGB>1)) = 1;
    figure(1)
    imagesc(RGB)
    axis equal
    
    LocalMaxThresh = InitialLocalMaxThresh;
    [Xfit1, Yfit1, Zfit1, Xgof1, Ygof1, Zgof1, Intensity1, Xwidth1, Ywidth1, Zwidth1] = fitMultipleFoci(ImageStack_left,LocalMaxThresh,300);
    NumLandMarks = length(Xfit1);
    while NumLandMarks<10 && LocalMaxThresh>15
        LocalMaxThresh = LocalMaxThresh-10;
        [Xfit1, Yfit1, Zfit1, Xgof1, Ygof1, Zgof1, Intensity1, Xwidth1, Ywidth1, Zwidth1] = fitMultipleFoci(ImageStack_left,LocalMaxThresh,300);
        NumLandMarks = length(Xfit1);
    end
    display([num2str(NumLandMarks) ' beads identified.']);
    % find the minimum distance to each bead (from all other beads): MinDisArray
    MinDisMatrix = ones(NumLandMarks,NumLandMarks)*ImageSize;
    for i = 1:NumLandMarks
        for j = 1:NumLandMarks
            if i~=j
                MinDisMatrix(i,j) = ((Xfit1(i)-Xfit1(j))^2+(Yfit1(i)-Yfit1(j))^2)^0.5;
            end
        end
    end
    MinDisArray = min(MinDisMatrix);
    
    LocalMaxThresh = InitialLocalMaxThresh;
    [Xfit2, Yfit2, Zfit2, Xgof2, Ygof2, Zgof2, Intensity2, Xwidth2, Ywidth2, Zwidth2] = fitMultipleFoci(ImageStack_right,LocalMaxThresh,300);
    NumLandMarks2 = length(Xfit2);
    while NumLandMarks2<10 && LocalMaxThresh>15
        LocalMaxThresh = LocalMaxThresh-10;
        [Xfit2, Yfit2, Zfit2, Xgof2, Ygof2, Zgof2, Intensity2, Xwidth2, Ywidth2, Zwidth2] = fitMultipleFoci(ImageStack_right,LocalMaxThresh,300);
        NumLandMarks2 = length(Xfit2);
    end
    display([num2str(NumLandMarks2) ' beads identified.']);
    % find the minimum distance to each bead (from all other beads): MinDisArray
    MinDisMatrix2 = ones(NumLandMarks2,NumLandMarks2)*ImageSize;
    for i = 1:NumLandMarks2
        for j = 1:NumLandMarks2
            if i~=j
                MinDisMatrix2(i,j) = ((Xfit2(i)-Xfit2(j))^2+(Yfit2(i)-Yfit2(j))^2)^0.5;
            end
        end
    end
    MinDisArray2 = min(MinDisMatrix2);
    
    % try to match the beads identified in later image to the beads
    % identified in the first image
    N = 0;
    for i = 1:NumLandMarks
        Xfit1_new = Xfit1(i);
        Yfit1_new = Yfit1(i);
        [m, Ind] = min(((Xfit2-Xfit1_new).^2+(Yfit2-Yfit1_new).^2).^0.5);
        if ~isempty(m)
            if m<MinDisArray(i)/2 && m<MinDisArray2(Ind)/2 && m<60
                N = N+1;
                Xfit1_match(N) = Xfit1(i);
                Yfit1_match(N) = Yfit1(i);
                Zfit1_match(N) = Zfit1(i);
                Xgof1_match(N) = Xgof1(i);
                Ygof1_match(N) = Ygof1(i);
                Zgof1_match(N) = Zgof1(i);
                Xfit2_match(N) = Xfit2(Ind);
                Yfit2_match(N) = Yfit2(Ind);
                Zfit2_match(N) = Zfit2(Ind);
                Xgof2_match(N) = Xgof2(Ind);
                Ygof2_match(N) = Ygof2(Ind);
                Zgof2_match(N) = Zgof2(Ind);
            end
        end
    end
    display([num2str(N) ' beads initially matched.']);
    
    
    MeanX = mean(Xfit1_match-Xfit2_match);
    MeanY = mean(Yfit1_match-Yfit2_match);
    A = [1, 0, 0; 0, 1, 0; MeanX, MeanY, 1];
    tform = affine2d(A);
    [x_new, y_new] = transformPointsForward(tform, Xfit2_match',Yfit2_match');
    MeanShiftX = median(Xfit1_match-x_new');
    MeanShiftY = median(Yfit1_match-y_new');
    Ind = find(abs(Xfit1_match-x_new'-MeanShiftX)<1 & abs(Yfit1_match-y_new'-MeanShiftY)<1);
    display([num2str(length(Ind)) ' beads finally matched.']);
    Xfit1_match = Xfit1_match(Ind);
    Yfit1_match = Yfit1_match(Ind);
    Zfit1_match = Zfit1_match(Ind);
    Xgof1_match = Xgof1_match(Ind);
    Ygof1_match = Ygof1_match(Ind);
    Zgof1_match = Zgof1_match(Ind);
    Xfit2_match = Xfit2_match(Ind);
    Yfit2_match = Yfit2_match(Ind);
    Zfit2_match = Zfit2_match(Ind);
    Xgof2_match = Xgof2_match(Ind);
    Ygof2_match = Ygof2_match(Ind);
    Zgof2_match = Zgof2_match(Ind);
    MeanX = mean(Xfit1_match-Xfit2_match);
    MeanY = mean(Yfit1_match-Yfit2_match);
    A = [1, 0, 0; 0, 1, 0; MeanX, MeanY, 1];
    tform = affine2d(A);
    [x_new, y_new] = transformPointsForward(tform, Xfit2_match',Yfit2_match');
    
    
    for i = 1:length(Xfit1_match)
        if Xgof1_match(i).adjrsquare>AdjrsquareThreshold && Ygof1_match(i).adjrsquare>AdjrsquareThreshold && Zgof1_match(i).adjrsquare>AdjrsquareThreshold ...
                && Xgof2_match(i).adjrsquare>AdjrsquareThreshold && Ygof2_match(i).adjrsquare>AdjrsquareThreshold && Zgof2_match(i).adjrsquare>AdjrsquareThreshold
            lmol_match = cat(1, lmol_match, [Xfit1_match(i), Yfit1_match(i), Zfit1_match(i)]);
            rmol_match = cat(1, rmol_match, [Xfit2_match(i), Yfit2_match(i), Zfit2_match(i)]);
        end
    end
    
end
%%
tform = cp2tform(rmol_match(:,1:2),lmol_match(:,1:2),'projective');
save('tform.mat','tform');
DeltaZ = mean(rmol_match(:,3)-lmol_match(:,3))
save('DeltaZ.mat','DeltaZ');
display([num2str(length(rmol_match)) ' beads finally matched in all FOVs.']);

%%
TransImg = imtransform(RGB(:,:,2), tform, 'XData', [1 ImageSize], 'Ydata', [1 ImageSize]);
RGB(:,:,2) = TransImg;
figure(2)
imagesc(RGB)
axis equal
%%
if exist('figures_calibration')
    delete figures_calibration/*
end
mkdir('figures_calibration')
figure(1)
savefig(['figures_calibration/fig1.fig']);
figure(2)
savefig(['figures_calibration/fig2.fig']);

