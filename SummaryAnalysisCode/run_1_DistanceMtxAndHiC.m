clear all
close all

TotalNumTADs = 28;
load('WT-IMR90\AllXaChr.mat')
load('WT-IMR90\AllXiChr.mat')
%%
% Fig 1B
Chr = [AllXaChr AllXiChr];

Mean_total = zeros(TotalNumTADs,TotalNumTADs);
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(Chr)
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
            end
        end
        Mean_total(i,j) = mean(DisList);
    end
end

figure
imagesc(Mean_total)
colorbar
caxis([0 0.7])
title(['Mean spatial distance, N=' num2str(k)]);
ColorMap = load('RedBlue.txt');
colormap(ColorMap/255);
PlotProp
axis square
% savefig('Figure1/Figure1_B.fig');
%%
% HiC plot
load('HiC.mat')

figure(8)
imagesc(HiC)
colorbar
caxis([0 650])
title('HiC interaction map');
axis square
ColorMap = load('RedBlue.txt');
L = size(ColorMap,1);
ColorMap = ColorMap(L:-1:1,:);
ColorMap = ColorMap(1:225,:)
N = 1000;
for i = 1:N
    ColorMap = [ColorMap; [192-(192-120)/N*i 0 0]];
end
ColorMap = ColorMap(1:6:end,:);
colormap(ColorMap/255);
PlotProp
%%
% HiC correlation
HiC_pxl = [];
Dis_pxl = [];
for i = 1:27
    for j = i+1:28
        HiC_pxl = [HiC_pxl HiC(i,j)];
        Dis_pxl = [Dis_pxl Mean_total(i,j)];
    end
end
Dis_pxl = Dis_pxl*1000;  % um to nm

figure(3)
loglog(Dis_pxl, HiC_pxl, '.','MarkerSize', 20);
set(gca, 'LineWidth', 1, 'FontSize', 16);
hold on
axis square
Ind = find(HiC_pxl>0 & Dis_pxl>0);
x = log(Dis_pxl(Ind))';
y = log(HiC_pxl(Ind))';
X = [ones(size(x)), x];
[b, bint] = regress(y,X);
text(280,10,['Slope = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))],'Color','red','FontSize',10)
R = corrcoef(x, y);
text(280,15,['Correlation coefficient = ' num2str(R(2,1))],'FontSize',10)
x = Dis_pxl';
y_fit = exp(b(1))*x.^b(2);
loglog(x, y_fit, '-r','LineWidth', 2);
hold off

%%
% plot Xa and Xi mean mtx
ChrXa = AllXaChr;
Mean_Xa = zeros(TotalNumTADs,TotalNumTADs);
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(ChrXa)
            if ChrXa(k).r(i) == 1 && ChrXa(k).r(j) == 1
                DisList = [DisList ((ChrXa(k).x(i)-ChrXa(k).x(j))^2+(ChrXa(k).y(i)-ChrXa(k).y(j))^2+(ChrXa(k).z(i)-ChrXa(k).z(j))^2)^0.5];
            end
        end
        Mean_Xa(i,j) = mean(DisList);
    end
end

ChrXi = AllXiChr;
Mean_Xi = zeros(TotalNumTADs,TotalNumTADs);
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(ChrXi)
            if ChrXi(k).r(i) == 1 && ChrXi(k).r(j) == 1
                DisList = [DisList ((ChrXi(k).x(i)-ChrXi(k).x(j))^2+(ChrXi(k).y(i)-ChrXi(k).y(j))^2+(ChrXi(k).z(i)-ChrXi(k).z(j))^2)^0.5];
            end
        end
        Mean_Xi(i,j) = mean(DisList);
    end
end

figure
imagesc(Mean_Xa)
colorbar
caxis([0 0.7])
title(['Xa Mean spatial distance, N=' num2str(length(ChrXa))]);
ColorMap = load('RedBlue.txt');
colormap(ColorMap/255);
PlotProp
axis square
figure
imagesc(Mean_Xi)
colorbar
caxis([0 0.7])
title(['Xi Mean spatial distance, N=' num2str(length(ChrXi))]);
ColorMap = load('RedBlue.txt');
colormap(ColorMap/255);
PlotProp
axis square
% savefig('Figure1/Figure1_E.fig');
%%
% figure 2, plot some examples of single-cell domain
clear all
close all
TotalNumTADs = 28;
load('WT-IMR90\AllXaChr.mat')
load('WT-IMR90\AllXiChr.mat')
% Xa
Chr = AllXaChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=23 % 80% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
Chr = Chr_chozen;
for k = 1:150
    Mean = ones(TotalNumTADs,TotalNumTADs);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                Mean(i,j) = ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5;
            end
        end
    end
    Idx = find(Mean==1);
    Mean(Idx) = -0.01;
    figure
    imagesc(Mean)
    colorbar
    caxis([-0.01 0.7])
    title(['Xa individual spatial distance-' num2str(k)]);
    ColorMap = load('RedBlue.txt');
    ColorMap = [[128 128 128]; ColorMap];
    colormap(ColorMap/255);
    PlotProp
    axis square
    savefig(['Figure2/Xa/Xa_' num2str(k) '.fig']);
end
close all

% Xi
Chr = AllXiChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=23 % 80% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
Chr = Chr_chozen;
for k = 1:150
    Mean = ones(TotalNumTADs,TotalNumTADs);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                Mean(i,j) = ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5;
            end
        end
    end
    Idx = find(Mean==1);
    Mean(Idx) = -0.01;
    figure
    imagesc(Mean)
    colorbar
    caxis([-0.01 0.7])
    title(['Xi individual spatial distance-' num2str(k)]);
    ColorMap = load('RedBlue.txt');
    ColorMap = [[128 128 128]; ColorMap];
    colormap(ColorMap/255);
    PlotProp
    axis square
    savefig(['Figure2/Xi/Xi_' num2str(k) '.fig']);
end
close all


















