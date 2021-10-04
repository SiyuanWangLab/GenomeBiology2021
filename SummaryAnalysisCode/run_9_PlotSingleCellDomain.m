clear all
close all

TotalNumTADs = 28;

SumFolder = 'SingleCellDomainFolder/';
mkdir(SumFolder);
Folder{1} = 'WT-IMR90';
Folder{2} = 'WT-TAD2';
Folder{3} = 'WT-TAD57';
Folder{4} = 'DMOG';
Folder{5} = 'GSK126';
Folder{6} = '5Aza';

Xchr{1} = 'Xi';
Xchr{2} = 'Xa';
n = 0;
for i = 1:length(Folder)
    if isempty(Folder{i})
        continue
    end
    for j=1:2
        n = n+1;
        FolderName{n} = [SumFolder Folder{i} '/' Xchr{j}];
        if ~exist(FolderName{n}, 'dir')
            mkdir(FolderName{n})
        end
    end
end
%%
for ii = 1:12
    
    Words = split(FolderName{ii},'/');
    ChrFile = [Words{2} '/All' Words{3} 'Chr'];
    display(ChrFile)
    load(ChrFile);
    Xchr = Words{3};
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
        if sum(Chr(k).r)>=23 % 80% or higher detection efficiency
            Chr_chozen = [Chr_chozen Chr(k)];
        end
    end
    Chr = Chr_chozen;
    
    for k = 1:50
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
        title([Xchr ' individual spatial distance-' num2str(k)]);
        ColorMap = load('RedBlue.txt');
        ColorMap = [[128 128 128]; ColorMap];
        colormap(ColorMap/255);
        PlotProp
        axis square
%         pause(1)
        savefig([FolderName{ii} '/Chr_' num2str(k) '.fig']);
        close
    end
    
end
