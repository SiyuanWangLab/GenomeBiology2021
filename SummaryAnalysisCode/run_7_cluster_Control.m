clear all
close all
TotalNumTADs = 28;
ChrName = 'WT';

load('WT-IMR90\AllXaChr.mat')
load('WT-IMR90\AllXiChr.mat')

rng default

Chr = AllXaChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=23 % 80% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
AllXaChr = Chr_chozen;

Chr = AllXiChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=23 % 80% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
AllXiChr = Chr_chozen;

Chr = [AllXaChr AllXiChr];
IdentityIdx = [ones(1,length(AllXaChr)) 2*ones(1,length(AllXiChr))];
%%
AllScChrDis = [];
RadiusOfGyration = [];
for k = 1:length(Chr)
    % plot the individual distance matrix of all chr
    Mean = zeros(TotalNumTADs,TotalNumTADs);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                Mean(i,j) = ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5;
            else
                Mean(i,j) = NaN;
            end
        end
    end
    Mean_filtered = fillmissing(Mean,'linear');
    Mean_filtered = fillmissing(Mean_filtered,'linear',2);
    
    ScChrDis = [];
    for row=1:TotalNumTADs-1
        for col=row+1:TotalNumTADs
            ScChrDis = [ScChrDis Mean_filtered(row,col)];
        end
    end
    AllScChrDis = [AllScChrDis;ScChrDis];
    
end
%%
LouvainCellTypeList = louvainJaccardClustering(AllScChrDis, 50); 
%%
% mean matrix for each cluster
for iii=1:max(LouvainCellTypeList)
    Idx = find(LouvainCellTypeList==iii);
    ClusterChr = Chr(Idx);
    
    Mean_total = zeros(TotalNumTADs,TotalNumTADs);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            DisList = [];
            for k = 1:length(ClusterChr)
                if ClusterChr(k).r(i) == 1 && ClusterChr(k).r(j) == 1
                    DisList = [DisList ((ClusterChr(k).x(i)-ClusterChr(k).x(j))^2+(ClusterChr(k).y(i)-ClusterChr(k).y(j))^2+(ClusterChr(k).z(i)-ClusterChr(k).z(j))^2)^0.5];
                end
            end
            Mean_total(i,j) = mean(DisList);
        end
    end
    
    figure
    imagesc(Mean_total)
    colorbar
    caxis([0 0.8])
    title(['Mean spatial distance, N=' num2str(k) ', cluster ' num2str(iii)]);
    ColorMap = load('RedBlue.txt');
    colormap(ColorMap/255);
    PlotProp
    axis square
    savefig(['ClusterAnalysis/MeanMtx/' num2str(iii) '.fig']);

end
%% Plot 2D tsne, umap and phate plots with the clustering results
% 2D TSNE
TSNEoutput = tsne(AllScChrDis,'Distance','cosine'); % the cosine distance functions as a normalization
%%
color = lines(8);

figure(1)
gscatter(TSNEoutput(:,1),TSNEoutput(:,2),LouvainCellTypeList,color,'.',12);
xlabel('tsne1');
ylabel('tsne2');
axis square

figure(2)
gscatter(TSNEoutput(:,1),TSNEoutput(:,2),IdentityIdx,'rb','.',12);
xlabel('tsne1');
ylabel('tsne2');
axis square
%%
% Xa=1, red; Xi=2, blue
ChrPct = [];
for iii=1:6
    Idx = find(LouvainCellTypeList==iii);
    ChrID = IdentityIdx(Idx);
    XaPct = length(find(ChrID==1))/length(ChrID);
    XiPct = length(find(ChrID==2))/length(ChrID);
    ChrPct = [ChrPct; [XaPct XiPct]];
end

figure(5)
ba = bar(ChrPct,'stacked','FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',1.5);
ba(1).CData = [1 0 0];
ba(2).CData = [0 0 1];
axis square


y = [50 50; 25 75; 30 70];
figure
ba = bar(y,'stacked', 'FaceColor','flat');
ba(1).CData = [0.3 0.3 0.7];
ba(2).CData = [1 1 1]*0.8;
































