% updated on 191106 to reduce the nearest neighbor threshold, and also
% calcualte the distance between new foci and the average foci position in
% a trace when trying to decided whether or not to add the new foci to the
% trace.

% updated on 190314.

clear all
close all
NFOV = 22; % number of fields of views
TotalTADNum = 28;
WidthThreshMax = 4; % this threshold removes fitted "foci" that are too wide.
WidthThreshMin = 1; % this threshold removes fitted "foci" that are too narrow.
AdjrsquareThreshold = 0.7;
ImageSize = 1536; % number of pxls
UmPerPxl = 0.108;
StepSize = 0.2; %um

NearestNeighborThreshold = 1; %um, distance threshold between new focus and the average position of all previous foci in the trace
TraceLengthThreshold = 15; % remove traces that are shorter than this length
%%
% the algorithm: Define each focus in hyb1 as the initiation point of one
% trace. For each focus in hyb1, find the closest focus in hyb2,
% check if the chozen hyb1 focus is the closest to the chozen hyb2 focus,
% if so, and if the two foci are close enough, link these two focus into
% one trace. If not, do not grow the trace, move on to the next focus in
% hyb1 (next trace). Finally, if there are remaining foci in hyb2, define
% them as the initiation points of additional traces.

if ~exist('Traces_SmallScale', 'dir')
    mkdir('Traces_SmallScale');
end
if ~exist('figs_Traces_SmallScale', 'dir')
    mkdir('figs_Traces_SmallScale');
end
FociCount = zeros(TotalTADNum,1); % initial foci count
FinalFociCount = zeros(TotalTADNum,1); % final count of foci in the traces

for jj = 0:NFOV-1
    
    if NFOV<=10
        FOVid = ['0' num2str(jj)];
    elseif NFOV>10 && NFOV<=100
        if jj<10
            FOVid = ['0' num2str(jj)];
        else
            FOVid = [num2str(jj)];
        end
    elseif NFOV>100
        if jj<10
            FOVid = ['00' num2str(jj)];
        elseif jj<100
            FOVid = ['0' num2str(jj)];
        else
            FOVid = [num2str(jj)];
        end
    end
    if exist(['TracingResults_SmallScale/result' FOVid '.mat'])==2
        load(['TracingResults_SmallScale/result' FOVid '.mat']);
    else
        continue
    end
   %% 
%     figure(1)
    for i = 1:TotalTADNum
        if ~isempty(XfitList{i})
            Ind = find(XwidthList{i}<WidthThreshMax & YwidthList{i}<WidthThreshMax ...
                & XwidthList{i}>WidthThreshMin & YwidthList{i}>WidthThreshMin ...
                & [XgofList{i}.adjrsquare]>AdjrsquareThreshold ...
                & [YgofList{i}.adjrsquare]>AdjrsquareThreshold ...
                & [ZgofList{i}.adjrsquare]>AdjrsquareThreshold);
            
            XfitList{i} = XfitList{i}(Ind);
            YfitList{i} = YfitList{i}(Ind);
            ZfitList{i} = ZfitList{i}(Ind);
            IntensityList{i} = IntensityList{i}(Ind);
            XwidthList{i} = XwidthList{i}(Ind);
            YwidthList{i} = YwidthList{i}(Ind);
            ZwidthList{i} = ZwidthList{i}(Ind);
            XgofList{i} = XgofList{i}(Ind);
            YgofList{i} = YgofList{i}(Ind);
            ZgofList{i} = ZgofList{i}(Ind);
            %                scatter(XfitList{i},YfitList{i},'.')
            %                hold on
        end
    end
    %         hold off
    %         axis equal
    %         axis ij
    for i = 1:TotalTADNum
        FociCount(i) = FociCount(i)+length(XfitList{i});
    end
    for i = 1:TotalTADNum
        P1 = [XfitList{i}'*UmPerPxl, YfitList{i}'*UmPerPxl, ZfitList{i}'*StepSize, ...
            IntensityList{i}', XwidthList{i}', YwidthList{i}', ZwidthList{i}', ...
            [XgofList{i}.adjrsquare]', [YgofList{i}.adjrsquare]', [ZgofList{i}.adjrsquare]'];
        if i == 1 % define initial traces
            for j = 1:size(P1,1)
                TraceArray{j} = [P1(j,:), i];
                %TraceArray: x(um),y(um),z(um), Intensity, Xwidth(pxl),
                %Ywidth(pxl), Zwidth(pxl), Xgof, Ygof, Zgof, hybNo
            end
        else % match to existing traces and define new traces
            P0 = [];
            for j = 1:length(TraceArray)
                % this following line was modified on 191106 to
                % calculate the average position of all foci in each
                % trace and add to P0
                P0 = [P0; mean(TraceArray{j}(:,1:3),1)]; % build the list of average foci positions in all traces
                %                     P0 = [P0; TraceArray{j}(end,1:3)]; % build the list of end foci in all traces
            end
            D = pdist2(P0, P1(:,1:3), 'euclidean');
            for j = 1:size(P0,1)
                [M, Ind] = min(D(j,:));
                M2 = min(D(:,Ind));
                if M<NearestNeighborThreshold && M==M2 && ~isnan(P1(Ind,1))
                    % add to existing trace
                    TraceArray{j}(end+1,:) = [P1(Ind,:),i];
                    % mark linked foci in P1
                    P1(Ind,1:3) = [NaN, NaN, NaN];
                end
            end
            % find remaining foci in P1 and define new traces
            for j = 1:size(P1,1)
                if ~isnan(P1(j,1))
                    TraceArray{end+1} = [P1(j,:), i]; % x,y,z,hybNo
                end
            end
        end
    end
    n = 0;
    for j = 1:length(TraceArray)
        if size(TraceArray{j},1) >= TraceLengthThreshold
            n = n+1;
            TraceArrayNew{n} = TraceArray{j};
        end
    end
    TraceArray = TraceArrayNew;
    clear TraceArrayNew
    
    figure(2)
    for j = 1:length(TraceArray)
        plot(TraceArray{j}(:,1), TraceArray{j}(:,2),'.-');
        hold on
        FinalFociCount(TraceArray{j}(:,end)) = FinalFociCount(TraceArray{j}(:,end))+1;
    end
    hold off
    axis equal
    axis ij
    xlabel('x (um)');
    ylabel('y (um)');
    savefig(['figs_Traces_SmallScale\Traces_' FOVid '.fig']);
    save(['Traces_SmallScale\TraceArray' FOVid '.mat'],'TraceArray');
    clear TraceArray
end
figure(3)
bar(FociCount)
xlabel('Hyb number')
ylabel('Initially identified foci count')
savefig(['figs_Traces_SmallScale/TracingFociCount_initial_SmallScale.fig']);
figure(4)
bar(FinalFociCount)
xlabel('Hyb number')
ylabel('Foci count after linking traces')
savefig(['figs_Traces_SmallScale/TracingFociCount_middle_SmallScale.fig']);