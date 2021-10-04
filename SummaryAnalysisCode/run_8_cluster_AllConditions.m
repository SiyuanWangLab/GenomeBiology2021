clear all
close all
TotalNumTADs = 28;

Folder{1} = 'WT-IMR90';
Folder{2} = 'TSA';
Folder{3} = 'amanitin';
Folder{4} = 'DMOG';
Folder{5} = 'GSK126';
Folder{6} = '5Aza';

rng default

ChrAll = [];
IdentityIdx = [];
N = 0;
for iii=1:length(Folder)
    load([Folder{iii} '\AllXaChr.mat'])
    load([Folder{iii} '\AllXiChr.mat'])
    Chr = [AllXaChr AllXiChr];
    Chr_chozen = [];
    for k = 1:length(Chr)
        if sum(Chr(k).r)>=23 % 80% or higher detection efficiency
            Chr_chozen = [Chr_chozen Chr(k)];
        end
    end
    Chr = Chr_chozen;
    for i=1:length(Chr)
        N=N+1;
        ChrAll(N).x = Chr(i).x;
        ChrAll(N).y = Chr(i).y;
        ChrAll(N).z = Chr(i).z;
        ChrAll(N).r = Chr(i).r;
    end
    n = length(Chr);
    Vec = iii*ones(1,n);
    IdentityIdx = [IdentityIdx Vec];
end

%%
AllScChrDis = [];
for k = 1:length(ChrAll)
    % plot the individual distance matrix of all chr
    Mean = zeros(TotalNumTADs,TotalNumTADs);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            if ChrAll(k).r(i) == 1 && ChrAll(k).r(j) == 1
                Mean(i,j) = ((ChrAll(k).x(i)-ChrAll(k).x(j))^2+(ChrAll(k).y(i)-ChrAll(k).y(j))^2+(ChrAll(k).z(i)-ChrAll(k).z(j))^2)^0.5;
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
LouvainCellTypeList = louvainJaccardClustering(AllScChrDis, 200); 
%%
% mean matrix for each cluster
for iii=1:max(LouvainCellTypeList)
    Idx = find(LouvainCellTypeList==iii);
    ClusterChr = ChrAll(Idx);
    
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
%     savefig(['ClusterAnalysis/Drug_MeanMtx/' num2str(iii) '.fig']);

end
%% Plot 2D tsne, umap and phate plots with the clustering results
% 2D TSNE
TSNEoutput = tsne(AllScChrDis,'Distance','cosine'); % the cosine distance functions as a normalization
%%
color = lines(6);
figure(1)
gscatter(TSNEoutput(:,1),TSNEoutput(:,2),LouvainCellTypeList,color,'.',12);
xlabel('tsne1');
ylabel('tsne2');
axis square

figure(2)
gscatter(TSNEoutput(:,1),TSNEoutput(:,2),IdentityIdx,color,'.',12);
xlabel('tsne1');
ylabel('tsne2');
axis square
%%
ChrPct = [];
% for iii=1:6 % condition
%     Idx = find(IdentityIdx==iii);
%     ClusterID = LouvainCellTypeList(Idx);
%     Vec = [];
%     for j=1:6 % 6 cluster
%         Pct = length(find(ClusterID==j))/length(Idx);
%         Vec = [Vec Pct];
%     end
%     ChrPct = [ChrPct; Vec];
% end
for iii=1:6
    Idx = find(LouvainCellTypeList==iii);
    ChrID = IdentityIdx(Idx);
    Vec = [];
    for j=1:6 % 6 condition
        Pct = length(find(ChrID==j))/length(find(IdentityIdx==j));
        Vec = [Vec Pct];
    end
    ChrPct = [ChrPct; Vec];
end
%%
figure(5)
ba = bar(ChrPct,'FaceColor','flat','EdgeColor',[0 0 0],'LineWidth',0.5);
axis square
legend()

































