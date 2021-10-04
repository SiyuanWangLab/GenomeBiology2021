clear all
close all
TotalNumTADs = 28;
load('WT-IMR90\AllXaChr.mat')
load('WT-IMR90\AllXiChr.mat')

DDD = 0.2; % Contact threshold = 200nm
%%
Chr = AllXiChr;
MatrixAll = [];
for k = 1:length(Chr)
    Mean = ones(28,28);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                Mean(i,j) = ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5;
            end
        end
    end
    MatrixAll{k} = Mean;
end

% T_Num = 0;
% Total_Num = 0;
P1 = [];
P2 = [];
P3 = [];
for i = 1:26
    for j = i+1:27
        for k = j+1:28
            
            ContactRecord = [];
            Num = 0;
            for n = 1:length(Chr)
                if sum([Chr(n).r(i) Chr(n).r(j) Chr(n).r(k)])<3
                    continue
                end
                Num = Num+1;
                CurrM = MatrixAll{n};
                ContactRecord{Num} = [CurrM(i,j)  CurrM(j,k) CurrM(i,k)];
            end
            C1 = cellfun(@(a) a(1,1), ContactRecord);
            C2 = cellfun(@(a) a(1,2), ContactRecord);
            C3 = cellfun(@(a) a(1,3), ContactRecord);
            
            p1 = length(find(C1<=DDD & C2<=DDD))/length(find(C1<=DDD));
            p2 = length(find(C2<=DDD))/length(C2);
            p3 = length(find(C1>DDD & C2<=DDD))/length(find(C1>DDD));
            P1 = [P1 p1];
            P2 = [P2 p2];
            P3 = [P3 p3];

        end
    end
end

[P2_Sort,Idx] = sort(P2);
P1_Sort = P1(Idx);
P3_Sort = P3(Idx);

X = 1:length(Idx);
figure(3)
scatter(X,P1_Sort,5,'filled','MarkerFaceColor','#ff3333')
hold on
scatter(X,P3_Sort,5,'filled','MarkerFaceColor','#4DBEEE')
hold on
plot(X,P2_Sort,'LineWidth',4,'Color','k');
% PlotProp
axis square
title('Xi')
%%
% Xa
Chr = AllXaChr;
MatrixAll = [];
for k = 1:length(Chr)
    Mean = ones(28,28);
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                Mean(i,j) = ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5;
            end
        end
    end
    MatrixAll{k} = Mean;
end

% T_Num = 0;
% Total_Num = 0;
P1 = [];
P2 = [];
P3 = [];
for i = 1:26
    for j = i+1:27
        for k = j+1:28
            
            ContactRecord = [];
            Num = 0;
            for n = 1:length(Chr)
                if sum([Chr(n).r(i) Chr(n).r(j) Chr(n).r(k)])<3
                    continue
                end
                Num = Num+1;
                CurrM = MatrixAll{n};
                ContactRecord{Num} = [CurrM(i,j)  CurrM(j,k) CurrM(i,k)];
            end
            C1 = cellfun(@(a) a(1,1), ContactRecord);
            C2 = cellfun(@(a) a(1,2), ContactRecord);
            C3 = cellfun(@(a) a(1,3), ContactRecord);
            
            p1 = length(find(C1<=DDD & C2<=DDD))/length(find(C1<=DDD));
            p2 = length(find(C2<=DDD))/length(C2);
            p3 = length(find(C1>DDD & C2<=DDD))/length(find(C1>DDD));
            P1 = [P1 p1];
            P2 = [P2 p2];
            P3 = [P3 p3];

        end
    end
end

[P2_Sort,Idx] = sort(P2);
P1_Sort = P1(Idx);
P3_Sort = P3(Idx);

X = 1:length(Idx);
figure(2)
scatter(X,P1_Sort,5,'filled','MarkerFaceColor','#ff3333')
hold on
scatter(X,P3_Sort,5,'filled','MarkerFaceColor','#4DBEEE')
hold on
plot(X,P2_Sort,'LineWidth',4,'Color','k');
% PlotProp
axis square
title('Xa')
















