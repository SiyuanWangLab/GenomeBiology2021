clear all
close all

TotalNumTADs = 28;

SumFolder = 'Figure4/Strength/';
mkdir(SumFolder);
Folder{1} = 'WT-IMR90';
Folder{2} = 'TSA';
Folder{3} = 'amanitin';
Folder{4} = 'DMOG';
Folder{5} = 'GSK126';
Folder{6} = '5Aza';

Xchr{1} = 'Xi';
Xchr{2} = 'Xa';
n = 0;
%%
for i = 1:length(Folder)
    for j=1:2
        n = n+1;
        FolderName{n} = [SumFolder Folder{i} '_' Xchr{j} '_Strength'];
        if ~exist(FolderName{n}, 'dir')
            mkdir(FolderName{n})
        end
    end
end
%%
for ii = 1:length(FolderName)
    
    Words = split(FolderName{ii},'_');
    Words_sub = split(Words{1},'/');
    ChrFile = [Words_sub{3} '/' 'All' Words{2} 'Chr'];
    load(ChrFile);
    Xchr = Words{2};
    clear Chr
    if strcmp(Xchr,'Xi')
        Chr = AllXiChr;
    elseif strcmp(Xchr,'Xa')
        Chr = AllXaChr;
    else
        display([ChrFile ' does not have Xi assignment!'])
    end
    Chr_chozen = [];
    for k = 1:length(Chr)
        if sum(Chr(k).r)>=23 % 90% or higher detection efficiency
            Chr_chozen = [Chr_chozen Chr(k)];
        end
    end
    Chr = Chr_chozen;
    
    for k = 1:length(Chr)
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
        
        
        for i = 1:TotalNumTADs
            if i >= 4 && i <= TotalNumTADs-8
                left(i) = median([Mean_filtered(i-3:i+6,i-3);Mean_filtered(i-2:i+7,i-2);Mean_filtered(i-1:i+8,i-1)]);
            elseif i > TotalNumTADs - 8
                left(i) = median([Mean_filtered(i-3:TotalNumTADs,i-3);Mean_filtered(i-2:TotalNumTADs,i-2);Mean_filtered(i-1:TotalNumTADs,i-1)]);
            elseif i == 1
                left(i) = Inf;
            elseif i == 2
                left(i) = median([Mean_filtered(i-1:TotalNumTADs,i-1)]);
            elseif i == 3
                left(i) = median([Mean_filtered(i-2:TotalNumTADs,i-2);Mean_filtered(i-1:TotalNumTADs,i-1)]);
            end
            
            if i <= TotalNumTADs - 11
                right(i) = median([Mean_filtered(i:i+9,i);Mean_filtered(i+1:i+10,i+1);Mean_filtered(i+2:i+11,i+2)]);
            elseif i > TotalNumTADs - 11 && i <= TotalNumTADs - 2
                right(i) = median([Mean_filtered(i:TotalNumTADs,i);Mean_filtered(i+1:TotalNumTADs,i+1);Mean_filtered(i+2:TotalNumTADs,i+2)]);
            elseif i == TotalNumTADs - 1
                right(i) = median([Mean_filtered(i:TotalNumTADs,i);Mean_filtered(i+1:TotalNumTADs,i+1)]);
            elseif i == TotalNumTADs
                right(i) = median([Mean_filtered(i:TotalNumTADs,i)]);
            end
            
            if i >= 12
                top(i) = median([Mean_filtered(i-11:i-2,i-2);Mean_filtered(i-10:i-1,i-1);Mean_filtered(i-9:i,i)]);
            elseif i >= 3 && i < 12
                top(i) = median([Mean_filtered(1:i-2,i-2);Mean_filtered(1:i-1,i-1);Mean_filtered(1:i,i)]);
            elseif i == 1
                top(i) = Mean_filtered(1:i,i);
            elseif i == 2
                top(i) = median([Mean_filtered(1:i-1,i-1);Mean_filtered(1:i,i)]);
            end
            
            if i >= 9 && i <= TotalNumTADs-3
                bottom(i) = median([Mean_filtered(i-6:i+3,i+3);Mean_filtered(i-8:i+1,i+1);Mean_filtered(i-7:i+2,i+2)]);
            elseif i < 9
                bottom(i) = median([Mean_filtered(1:i+3,i+3);Mean_filtered(1:i+1,i+1);Mean_filtered(1:i+2,i+2)]);
            elseif i == TotalNumTADs - 1
                bottom(i) = median(Mean_filtered(i-8:i+1,i+1));
            elseif i == TotalNumTADs
                bottom(i) = Inf;
            elseif i == TotalNumTADs - 2
                bottom(i) = median([Mean_filtered(1:i+1,i+1);Mean_filtered(1:i+2,i+2)]);
            end
            start(i) = left(i)/right(i);
            ending(i) = bottom(i)/top(i);
        end
        save([FolderName{ii} '/StartAndEnd_Chr_' num2str(k) '.mat'],'start','ending');
    end
    
    
end 
