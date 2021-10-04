function [ImageStack, InfoFile] = ReadZStack_MultiChannel_Trigger(FileName, NumImage, StepInterval,TotalNumChannels, ChannelIndex)
% 201023: This function was derived from ReadZStack_MultiChannel after adding the
% trigger cable to connect the camera to the DAQ board for both Jellyfish
% (201023) and Nautilus (201022). Both systems now have the same timing.
% For a three channel z stack with 15 frames at each height and 5 frames for
% each channel, we can averaged frames 5-7, 9-12, 14-16 respectively for 
% 647, 560 and 488 channels. Then the z positioner movement happens on frame 17.

[MovieFP, InfoFile] = ReadDax([FileName, '.dax'],'startFrame', 1, 'endFrame', NumImage);

% update on 2021.07.19, need to change the frame to match each channel,
% average only 3 frames
NumHeights = floor(NumImage/(StepInterval*TotalNumChannels)); % number of heights in z stack
StartingFrames = (((1:NumHeights)-1)*TotalNumChannels+ChannelIndex-1)*StepInterval+4; % starting frames for each averaging
EndingFrames = (((1:NumHeights)-1)*TotalNumChannels+ChannelIndex)*StepInterval+1; % ending frames for each averaging
if ChannelIndex == TotalNumChannels
    EndingFrames = EndingFrames-1; % for the last channel, or for single-channel imaging
elseif ChannelIndex == 1
    StartingFrames = StartingFrames +1; % for the first channel in multi-channel imaging
end

for i = 1:NumHeights
    ImageStack(:,:,i) = mean(MovieFP(:,:, StartingFrames(i):EndingFrames(i)),3);
end


