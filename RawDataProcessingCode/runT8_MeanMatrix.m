clear all
close all
TotalNumTADs = 28;
DomainsToExclude = []; % list the domains to exclude, e.g. DomainsToExclude = [15, 2];

load('AllXiChr.mat');
load('AllXaChr.mat');
load('AllXXChr.mat');

%%
for ii = 1:3
    if ii == 1
        Chr = AllXiChr;
        ChrName = 'Xi';
    elseif ii == 2
        Chr = AllXaChr;
        ChrName = 'Xa';
    else
        Chr = AllXXChr;
        ChrName = 'XX';
    end
    Mean = zeros(TotalNumTADs,TotalNumTADs);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            DisList = [];
            for k = 1:length(Chr)
                if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                    DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
                end
            end
            Mean(i,j) = mean(DisList);
%             Std(i,j) = std(DisList);
%             SEM(i,j) = std(DisList)/(length(DisList))^0.5;
%             NofData(i,j) = length(DisList);
%             DisListAll{i,j} = DisList;
        end
    end
    
    figure(ii)
    imagesc(Mean)
    colorbar
    title([ChrName ', Mean spatial distance, ' num2str(length(Chr)) ' chr']);
    caxis([0 0.8])
    ColorMap = load('RedBlue.txt');
    colormap(ColorMap/255);
    PlotProp
    axis square
    savefig([ChrName 'Mean spatial distance all chr.fig']);
end





















