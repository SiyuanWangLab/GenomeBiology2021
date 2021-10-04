clear all
close all

TotalNumTADs = 28;

SumFolder = 'Figure4/Strength/';
% mkdir(SumFolder);
Folder{1} = 'WT-IMR90';
Folder{2} = 'TSA';
Folder{3} = 'DMOG';
Folder{4} = 'GSK126';
Folder{5} = 'amanitin';
Folder{6} = '5Aza';

Xchr{1} = 'Xi';
Xchr{2} = 'Xa';
n = 0;
%%
for i = 1:length(Folder)
    for j=1:2
        n = n+1;
        FolderName{n} = [SumFolder Folder{i} '_' Xchr{j} '_Strength'];
    end
end
%%
ChrNum = [];
StartNum = [];
EndNum = [];
for ii = 1:length(FolderName)
    
    Words = split(FolderName{ii},'_');
    Words_sub = split(Words{1},'/');
    ChrFile = [Words_sub{3} '/' 'All' Words{2} 'Chr'];
    load(ChrFile);
    Xchr = Words{2};
    CurrFolder = FolderName{ii};
    
    clear start_pks start_locs ending_pks ending_locs BoundaryID
    
    AllFiles = dir([CurrFolder '/*.mat']);
    StartStrength = [];
    EndStrength = [];
    for j = 1:length(AllFiles)
        load([CurrFolder '/' AllFiles(j).name]);
        start = start(3:TotalNumTADs-2);
        ending = ending(3:TotalNumTADs-2);
        [start_pks{j},start_locs{j}] = findpeaks(start,'MinPeakHeight',1.8);
        [ending_pks{j},ending_locs{j}] = findpeaks(ending,'MinPeakHeight',1.8);
        StartStrength = [StartStrength start_pks{j}];
        EndStrength = [EndStrength ending_pks{j}];
        TADID = sort([start_locs{j} ending_locs{j}]);
%         Idx = [];
%         for i=1:length(TADID)
%             if isempty(Idx)
%                 Idx = [Idx TADID(i)];
%             elseif (TADID(i)-Idx(end))>=2
%                 Idx = [Idx TADID(i)];
%             end
%         end
%         BoundaryID{j} = Idx;
        BoundaryID{j} = TADID;
    end
    StrengthAll{ii} = [StartStrength EndStrength];
    Idx = find(~isinf(StrengthAll{ii}));
    StrengthAll{ii} = StrengthAll{ii}(Idx);
    BoundaryNumAll{ii} = cellfun(@length,BoundaryID)/2;
    BoundaryNumTotal(ii) = length(cell2mat(BoundaryID))/2;
    Condition{ii} = [Words_sub{3} ' ' Words{2}];
    ChrNum(ii)=length(AllFiles);
    
end
%%
% for i = 1:length(BoundaryNumAll)
%     ErrorBar(i) = std(BoundaryNumAll{i})/sqrt(length(BoundaryNumAll{i}));
% end
for i = 1:length(BoundaryNumAll)
    ErrorBar(i) = std(BoundaryNumAll{i});
end
BoundaryFrequencyAll = BoundaryNumTotal./ChrNum;
Num = length(Condition);
Boundary_Organize = [BoundaryFrequencyAll(1:2:Num); BoundaryFrequencyAll(2:2:Num)]';
ErrorBar = [ErrorBar(1:2:Num); ErrorBar(2:2:Num)]';

%%
Names = Folder;
figure
b = bar(Boundary_Organize);
hold on
% b(1).CData(5,:) = [0 1 0]; % group 1 1st bar
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(Boundary_Organize);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',Boundary_Organize,ErrorBar,'k','linestyle','none');
set(gca,'xticklabel',Names)
legend('Xi','Xa','Location','northeast');
hold off

%%
% Xi t-test
Xi_P = [];
[h,p] = ttest2(BoundaryNumAll{1},BoundaryNumAll{3});
Xi_P = [Xi_P p];
[h,p] = ttest2(BoundaryNumAll{1},BoundaryNumAll{5});
Xi_P = [Xi_P p];
[h,p] = ttest2(BoundaryNumAll{1},BoundaryNumAll{7});
Xi_P = [Xi_P p];
[h,p] = ttest2(BoundaryNumAll{1},BoundaryNumAll{9});
Xi_P = [Xi_P p];
[h,p] = ttest2(BoundaryNumAll{1},BoundaryNumAll{11});
Xi_P = [Xi_P p]

% Xa t-test
Xa_P = [];
[h,p] = ttest2(BoundaryNumAll{2},BoundaryNumAll{4});
Xa_P = [Xa_P p];
[h,p] = ttest2(BoundaryNumAll{2},BoundaryNumAll{6});
Xa_P = [Xa_P p];
[h,p] = ttest2(BoundaryNumAll{2},BoundaryNumAll{8});
Xa_P = [Xa_P p];
[h,p] = ttest2(BoundaryNumAll{2},BoundaryNumAll{10});
Xa_P = [Xa_P p];
[h,p] = ttest2(BoundaryNumAll{2},BoundaryNumAll{12});
Xa_P = [Xa_P p]

p_adjust_Xi = mafdr(Xi_P,'BHFDR', true)
p_adjust_Xa = mafdr(Xa_P,'BHFDR', true)
%%
% plot figure F
AllLength = cellfun(@length,StrengthAll);
MaxLength = max(AllLength);
StrengthXi = NaN(MaxLength,ngroups);
StrengthXa = NaN(MaxLength,ngroups);
for i = 1:ngroups
    StrengthXi(1:AllLength(2*i-1),i) = StrengthAll{2*i-1};
    StrengthXa(1:AllLength(2*i),i) = StrengthAll{2*i};
end

X1=1:ngroups;
X2=X1+0.3;
figure
boxplot(StrengthXi(:,1:ngroups),'positions',X1,'labels',X1,'Notch','on','colors','b','symbol', '','width',0.25)
ylim([1.5 4.5])
xlim([0.5 7])
hold on
boxplot(StrengthXa(:,1:ngroups),'positions',X2,'labels',X2,'Notch','on','colors','r','symbol', '','width',0.25)
ylim([1.5 4.7])
xlim([0.5 7])
hold on
set(gca,'xticklabel',Folder)
title('Boundary Strength Distribution')

%%
% Xi t-test
Xi_P = [];
[h,p] = ttest2(StrengthAll{1},StrengthAll{3});
Xi_P = [Xi_P p];
[h,p] = ttest2(StrengthAll{1},StrengthAll{5});
Xi_P = [Xi_P p];
[h,p] = ttest2(StrengthAll{1},StrengthAll{7});
Xi_P = [Xi_P p];
[h,p] = ttest2(StrengthAll{1},StrengthAll{9});
Xi_P = [Xi_P p];
[h,p] = ttest2(StrengthAll{1},StrengthAll{11});
Xi_P = [Xi_P p];


% Xa t-test
Xa_P = [];
[h,p] = ttest2(StrengthAll{2},StrengthAll{4});
Xa_P = [Xa_P p];
[h,p] = ttest2(StrengthAll{2},StrengthAll{6});
Xa_P = [Xa_P p];
[h,p] = ttest2(StrengthAll{2},StrengthAll{8});
Xa_P = [Xa_P p];
[h,p] = ttest2(StrengthAll{2},StrengthAll{10});
Xa_P = [Xa_P p];
[h,p] = ttest2(StrengthAll{2},StrengthAll{12});
Xa_P = [Xa_P p];

p_adjust_Xi = mafdr(Xi_P,'BHFDR', true)
p_adjust_Xa = mafdr(Xa_P,'BHFDR', true)











