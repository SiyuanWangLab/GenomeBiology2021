clear all
close all
TotalNumTADs = 28;
ChrName = 'ChrX';
NFOV = 22; % number of fields of views

%% analyze mean spatial distance matrix
n = 0;
for jj = 1:NFOV-1
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
    if exist(['Traces_SmallScale\TraceArrayRefined' FOVid '.mat'])==2
        load(['Traces_SmallScale\TraceArrayRefined' FOVid '.mat']);
        for i = 1:length(TraceArray)
            n = n+1;
            Chr(n).x = zeros(TotalNumTADs,1);
            Chr(n).y = zeros(TotalNumTADs,1);
            Chr(n).z = zeros(TotalNumTADs,1);
            Chr(n).r = zeros(TotalNumTADs,1);
            Chr(n).x(TraceArray{i}(:,end)) = TraceArray{i}(:,1);
            Chr(n).y(TraceArray{i}(:,end)) = TraceArray{i}(:,2);
            Chr(n).z(TraceArray{i}(:,end)) = TraceArray{i}(:,3);
            Chr(n).r(TraceArray{i}(:,end)) = 1;
        end
    end
end
display([num2str(n) ' copies of Chr' ChrName ' were analyzed.'])
%%
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(Chr)
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
            end
        end
        Mean(i,j) = mean(DisList);
%         Std(i,j) = std(DisList);
%         SEM(i,j) = std(DisList)/(length(DisList))^0.5;
%         NofData(i,j) = length(DisList);
%         DisListAll{i,j} = DisList;
    end
end

figure(2)
imagesc(Mean)
colorbar
title(['All Chr Mean spatial distance, ' num2str(length(Chr)) ' chr']);
caxis([0 0.8])
ColorMap = load('RedBlue.txt');
colormap(ColorMap/255);
PlotProp
axis square
savefig(['Mean spatial distance all chr.fig']);





















