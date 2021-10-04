clear all
close all

TotalNumTADs = 28;
load('WT-IMR90\AllXaChr.mat')
load('WT-IMR90\AllXiChr.mat')
%%
% Xi
Chr = AllXiChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=25 % 80% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
Chr = Chr_chozen;
% for each chromatin, calculate radius of gyration
RadiusOfGyration = [];
for k = 1:length(Chr)
    CurrChr = Chr(k);
    NotNaIdx = find(CurrChr.r==1);
    X = CurrChr.x(NotNaIdx);
    Y = CurrChr.y(NotNaIdx);
    Z = CurrChr.z(NotNaIdx);
    XCenter = mean(X);
    YCenter = mean(Y);
    ZCenter = mean(Z);
    R = [];
    for j = 1:length(NotNaIdx)
        R(j) = (X(j)-XCenter)^2+(Y(j)-YCenter)^2+(Z(j)-ZCenter)^2;
    end
    Rg = (sum(R)/length(R))^0.5;
    RadiusOfGyration = [RadiusOfGyration Rg];
end
RadiusOfGyrationAll{1} = RadiusOfGyration;

% Xa
Chr = AllXaChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=25 % 90% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
Chr = Chr_chozen;
% for each chromatin, calculate radius of gyration
RadiusOfGyration = [];
for k = 1:length(Chr)
    CurrChr = Chr(k);
    NotNaIdx = find(CurrChr.r==1);
    X = CurrChr.x(NotNaIdx);
    Y = CurrChr.y(NotNaIdx);
    Z = CurrChr.z(NotNaIdx);
    XCenter = mean(X);
    YCenter = mean(Y);
    ZCenter = mean(Z);
    R = [];
    for j = 1:length(NotNaIdx)
        R(j) = (X(j)-XCenter)^2+(Y(j)-YCenter)^2+(Z(j)-ZCenter)^2;
    end
    Rg = (sum(R)/length(R))^0.5;
    RadiusOfGyration = [RadiusOfGyration Rg];
end
RadiusOfGyrationAll{2} = RadiusOfGyration;
%%
AllLength = cellfun(@length,RadiusOfGyrationAll)
MaxLength = max(AllLength);
RGAll = NaN(MaxLength,2);
RGAll(1:AllLength(1),1) = RadiusOfGyrationAll{1};
RGAll(1:AllLength(2),2) = RadiusOfGyrationAll{2};
%%
figure(1)
boxplot(RGAll,'Notch','on','colors','br','symbol','','width',0.35)
ylim([0.1 0.7])
set(gca,'xticklabel',{'Xi','Xa'})
axis square
title('Small Scale Radius Of Gyration')

% [h,p] = ttest2(RadiusOfGyrationAll{1},RadiusOfGyrationAll{2})
%%
% clear all
% close all
TotalNumTADs = 40;
load('StevenData\AllXaChr.mat')
load('StevenData\AllXiChr.mat')
% Xi
Chr = AllXiChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=35 % 90% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
Chr = Chr_chozen;
% for each chromatin, calculate radius of gyration
RadiusOfGyration = [];
for k = 1:length(Chr)
    CurrChr = Chr(k);
    NotNaIdx = find(CurrChr.r==1);
    X = CurrChr.x(NotNaIdx);
    Y = CurrChr.y(NotNaIdx);
    Z = CurrChr.z(NotNaIdx);
    XCenter = mean(X);
    YCenter = mean(Y);
    ZCenter = mean(Z);
    R = [];
    for j = 1:length(NotNaIdx)
        R(j) = (X(j)-XCenter)^2+(Y(j)-YCenter)^2+(Z(j)-ZCenter)^2;
    end
    Rg = (sum(R)/length(R))^0.5;
    RadiusOfGyration = [RadiusOfGyration Rg];
end
RadiusOfGyrationAll{1} = RadiusOfGyration;

% Xa
Chr = AllXaChr;
Chr_chozen = [];
for k = 1:length(Chr)
    if sum(Chr(k).r)>=35 % 90% or higher detection efficiency
        Chr_chozen = [Chr_chozen Chr(k)];
    end
end
Chr = Chr_chozen;
% for each chromatin, calculate radius of gyration
RadiusOfGyration = [];
for k = 1:length(Chr)
    CurrChr = Chr(k);
    NotNaIdx = find(CurrChr.r==1);
    X = CurrChr.x(NotNaIdx);
    Y = CurrChr.y(NotNaIdx);
    Z = CurrChr.z(NotNaIdx);
    XCenter = mean(X);
    YCenter = mean(Y);
    ZCenter = mean(Z);
    R = [];
    for j = 1:length(NotNaIdx)
        R(j) = (X(j)-XCenter)^2+(Y(j)-YCenter)^2+(Z(j)-ZCenter)^2;
    end
    Rg = (sum(R)/length(R))^0.5;
    RadiusOfGyration = [RadiusOfGyration Rg];
end
RadiusOfGyrationAll{2} = RadiusOfGyration;
%%
AllLength = cellfun(@length,RadiusOfGyrationAll)
MaxLength = max(AllLength);
RGAll = NaN(MaxLength,2);
RGAll(1:AllLength(1),1) = RadiusOfGyrationAll{1};
RGAll(1:AllLength(2),2) = RadiusOfGyrationAll{2};
%%
figure(2)
boxplot(RGAll,'Notch','on','colors','rb','symbol','','width',0.35)
axis square
set(gca,'xticklabel',{'Xi','Xa'})
title('Large Scale Radius Of Gyration')

% [h,p] = ttest2(RadiusOfGyrationAll{1},RadiusOfGyrationAll{2})

























