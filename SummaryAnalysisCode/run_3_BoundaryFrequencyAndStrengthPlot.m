clear all
close all

TotalNumTADs = 28;

%%
% Xa
ChrNum = [];
StartNum = [];
EndNum = [];
CurrFolder = 'Figure2/IMR90/StrengthXa';
clear start_pks start_locs ending_pks ending_locs
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
%     StartStrength(start_locs{j}) = StartStrength(start_locs{j})+start_pks{j};
%     EndStrength(ending_locs{j}) = EndStrength(ending_locs{j})+ending_pks{j};
end
N_Xa = length(start_pks)
StartAll = cell2mat(start_locs);
EndAll = cell2mat(ending_locs);
Num_start = zeros(1,TotalNumTADs-4);
Num_end = zeros(1,TotalNumTADs-4);
for i=1:TotalNumTADs-4
    Num = length(find(StartAll==i));
    Num_start(i) = Num;
    Num = length(find(EndAll==i));
    Num_end(i) = Num;
end
for i=1:TotalNumTADs-4
    Probability_start(i) = Num_start(i)/length(start_locs)*100;
    Probability_end(i) = Num_end(i)/length(start_locs)*100;
    Domain_probability(i) = (Probability_start(i) + Probability_end(i))/2;
end

Line1 = mean(Domain_probability);
Line2 = Line1-std(Domain_probability);
Line3 = Line1+std(Domain_probability);
DP_Xa = Domain_probability;

figure(1)
plot(3:26,Domain_probability,'black','LineWidth',1,'Marker','o','MarkerEdgeColor','black','MarkerFaceColor',[0.5,0.5,0.5]);
hold on
ylim([0 12])
hline1 = refline([0 Line1]);
hline2 = refline([0 Line2]);
hline3 = refline([0 Line3]);
set(gca,'box','off','FontSize',15);
xlabel('Domain ID','FontSize',15);
title('Xa Boundary Probability','FontSize',10);
% savefig('Figure2/IMR90/Fig2D.fig');

StrengthAll = [StartStrength EndStrength];
Idx = find(~isinf(StrengthAll));
StrengthAll = StrengthAll(Idx);

BinSize = 20;
Density = zeros(1,BinSize);
X = [];
for i = 1:BinSize
    Start = 1+0.2*(i-1);
    End = 1+0.2*i;
    X = [X (Start+End)/2];
    Idx = find(StrengthAll>=Start & StrengthAll<End);
    Density(i) = length(Idx)/length(StrengthAll);
end
figure(2)
B = bar(X,Density,'histc');
B.FaceColor = '#05C603';
B.LineWidth = 1;
ylim([0 0.35])


%%
% Xi
ChrNum = [];
StartNum = [];
EndNum = [];
CurrFolder = 'Figure2/IMR90/StrengthXi';
clear start_pks start_locs ending_pks ending_locs
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
%     StartStrength(start_locs{j}) = StartStrength(start_locs{j})+start_pks{j};
%     EndStrength(ending_locs{j}) = EndStrength(ending_locs{j})+ending_pks{j};
end
N_Xi = length(start_pks)
StartAll = cell2mat(start_locs);
EndAll = cell2mat(ending_locs);
Num_start = zeros(1,TotalNumTADs-4);
Num_end = zeros(1,TotalNumTADs-4);
for i=1:TotalNumTADs-4
    Num = length(find(StartAll==i));
    Num_start(i) = Num;
    Num = length(find(EndAll==i));
    Num_end(i) = Num;
end
for i=1:TotalNumTADs-4
    Probability_start(i) = Num_start(i)/length(start_locs)*100;
    Probability_end(i) = Num_end(i)/length(start_locs)*100;
    Domain_probability(i) = (Probability_start(i) + Probability_end(i))/2;
end

Line1 = mean(Domain_probability);
Line2 = Line1-std(Domain_probability);
Line3 = Line1+std(Domain_probability);
DP_Xi = Domain_probability;

figure(11)
plot(3:26,Domain_probability,'black','LineWidth',1,'Marker','o','MarkerEdgeColor','black','MarkerFaceColor',[0.5,0.5,0.5]);
hold on
ylim([0 12])
hline1 = refline([0 Line1]);
hline2 = refline([0 Line2]);
hline3 = refline([0 Line3]);
set(gca,'box','off','FontSize',15);
xlabel('Domain ID','FontSize',15);
title('Xi Boundary Probability','FontSize',10);
% savefig('Figure2/IMR90/Fig2C.fig');


StrengthAll = [StartStrength EndStrength];
Idx = find(~isinf(StrengthAll));
StrengthAll = StrengthAll(Idx);

BinSize = 20;
Density = zeros(1,BinSize);
X = [];
for i = 1:BinSize
    Start = 1+0.2*(i-1);
    End = 1+0.2*i;
    X = [X (Start+End)/2];
    Idx = find(StrengthAll>=Start & StrengthAll<End);
    Density(i) = length(Idx)/length(StrengthAll);
end
figure(6)
B = bar(X,Density,'histc');
B.FaceColor = '#D93703';
B.LineWidth = 1;
ylim([0 0.35])
%%
% F-test
% DP_Xi
% DP_Xa
[h,p] = vartest2(DP_Xi,DP_Xa)
N_Xa
N_Xi



