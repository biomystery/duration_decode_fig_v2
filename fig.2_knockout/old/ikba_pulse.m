% Pull in and make an attempt to clean up the pulse experiments (pulsed LPS in aKO cells was particularly messy due to
% poor nuclear staining).
dir1 = '/Volumes/data-1/backed_up/brooks/Code/cell tracking/';
addpath('/Volumes/data-1/backed_up/brooks/MACKtrack/Utilities/')
addpath(genpath('/Volumes/data-1/backed_up/brooks/Code/information theory'))
setcolors;
savedir = './';
if ~exist(savedir,'dir'); mkdir(savedir); end

%colors.lps = colors.blue;
clr = {colors.lps, colors.lps, [0.2667    0.7353    0.4667],[0.2667    0.7353    0.4667]};


if ~exist([dir1, 'ikbapulse.mat'],'file')
    ids = 365:368;
    names = {'wt LPS', 'aKO LPS', 'aKO TNF', 'wt TNF'};
    pulse = struct;
    for i = 1:length(ids)
        pulse(i).metrics = nfkbmetrics(ids(i));
        pulse(i).ids = ids;
        pulse(i).name = names{i};
        pulse(i).id = ids(i);
    end
    save([dir1, 'ikbapulse.mat'],'pulse')
    clear i ids
else
    load([dir1, 'ikbapulse.mat'])
end
%%
mod_colormap = divergingmap(0:1/1023:1,[32 19 115]/255,[138 47 20]/255);
mod_colormap(1,:) = [0.1 0.1 0.1];

clr1 = [34 50 100];
cmap = [linspace(255,clr1(1),128)',linspace(255,clr1(2),128)',linspace(255,clr1(3),128)']/255;
%cmap = [0.2*ones(1,3);cmap];

% For each condition, show single cell traces
t = 0:(1/12):8;
for i = 1:length(pulse(1).ids)
    drops = nansum(pulse(i).metrics.time_series<0.2,2)>50;
    nfkb = smoothrows(pulse(i).metrics.time_series(~drops,:),5);
    [~,sortorder] = sort(nansum(nfkb(:,1:60),2),'descend');
    if i==1; 
        nfkb=nfkb(sortorder,1:97);
    else
        nfkb = nfkb(sortorder,2:98);
    end
    nfkb = nfkb(round(linspace(1,size(nfkb,1),350)),:);
    figs.heatmap(i) = figure('PaperPositionMode','auto', 'Position', [40 40 200 400],'Name', pulse(i).name);
    colormap(cmap)
    imagesc([min(t), max(t)], [1 size(nfkb,1)], nfkb,[0.3 4.5])
    set(gca,'XTickLabel',{},'YTickLabel',{},'XTick',[0 3 6])
end

figs.line_plot = figure('PaperPositionMode','auto', 'Position', ...
    [40 40 600 400],'Name','Population averages');
hold on
linestyles = {'-', '--', '--', '-'};
for i = 1:length(pulse(1).ids)
    drops = nansum(pulse(i).metrics.time_series<0.2,2)>25;
    nfkb = smoothrows(pulse(i).metrics.time_series(~drops,:),5);
    [~,sortorder] = sort(nansum(nfkb(:,1:60),2),'descend');
    if i==1; 
        nfkb=nfkb(sortorder,1:97);
    else
        nfkb = nfkb(sortorder,2:98);
    end
    plot(t, nanmean(nfkb),'LineWidth',2,'Color',clr{i},'Linestyle',linestyles{i})
end
hold off

%% show number of cells showing in each figure 
t = 0:(1/12):8;
for i = 1:length(pulse(1).ids)
    disp(pulse(i).name)
    disp('total of cells') 
    disp(size(pulse(i).metrics.time_series,1))
    disp('Dropped cells:')
    drops = nansum(pulse(i).metrics.time_series<0.2,2)>50;
    disp(sum(drops))
    disp('Shown cells:')
    disp(size(pulse(i).metrics.time_series(~drops,:),1))
    
end


%%
figs.colorbar = figure('Position',...
    [40 40 30 500],'paperpositionmode','auto');
imagesc((64:-1:1)'),colormap(cmap)
set(gca,'Xtick',[],'YTick',[])
%% Mututal information plots

for tp = 16%13:1:20
    timepts = [8 tp];

    disp(['Timepoints: [',num2str((timepts(1)-1)*5),',',num2str((timepts(2)-1)*5),']'])
    figs.state_space = figure('PaperPositionMode','auto', 'Position', ...
        [40 40 650 280],'Name','Population averages');
    ha = tight_subplot(1,2,[0.01 0.07],[0.08 0.05]);
    % WT cells:
    hold(ha(1),'on')
    X= cell(2,1);
    idx = [1 4];
    for i = 1:length(idx)
        drops = nansum(pulse(idx(i)).metrics.time_series<0.2,2)>50;
        nfkb = smoothrows(pulse(idx(i)).metrics.time_series(~drops,1:length(t)),5);
        if i==1; 
            nfkb=nfkb(:,1:96);
        else
            nfkb = nfkb(:,2:97);
        end
        X{i} = nfkb(:,timepts)';
        % Downsample NFkB so similar numbers of cells are plotted
        nfkb = nfkb(1:min([350,size(nfkb,1)]),:);
        plot(ha(1),nfkb(:,timepts(1)), nfkb(:,timepts(2)),'.','MarkerSize',12, 'Color', clr{idx(i)})%[clr{idx(i)} 0.4])
        set(ha(1),'xlim',[0 12], 'ylim',[0 12],'xtick',0:3:12, 'ytick',0:3:12,'box','on')
    end
    [fitI,I,~,Q] = jacknifeCC(X,9,12);       
    disp(['wt I=', num2str(fitI)])
    text(0.3,11.7,  ['MI = ', num2str(round(fitI*100)/100),' bits'],'Parent',ha(1),'FontSize',14,'VerticalAlignment','top')
    
    % aKO cells:
    hold(ha(2),'on')
    X= cell(2,1);
    idx = [2 3];
    for i = 1:length(idx)
        drops = nansum(pulse(idx(i)).metrics.time_series<0.2,2)>50;
        nfkb = smoothrows(pulse(idx(i)).metrics.time_series(~drops,1:length(t)),5);
        if i==1; 
            nfkb=nfkb(:,1:96);
        else
            nfkb = nfkb(:,2:97);
        end
        X{i} = nfkb(:,timepts)';
        % Downsample NFkB so similar numbers of cells are plotted
        nfkb = nfkb(1:min([350,size(nfkb,1)]),:);
        plot(ha(2),nfkb(:,timepts(1)), nfkb(:,timepts(2)),'.','MarkerSize',12, 'Color', clr{idx(i)})%[clr{idx(i)} 0.4])
        set(ha(2),'xlim',[0 8], 'ylim',[0 8],'xtick',0:2:8, 'ytick',0:2:8,'box','on')
        % Add in date for information run
    end
    [fitI,I,~,Q] = jacknifeCC(X,9,12);  
    text(0.2, 7.8, ['MI = ', num2str(round(fitI*100)/100),' bits'],'Parent',ha(2),'FontSize',14,'VerticalAlignment','top')

    disp(['aKO I=', num2str(fitI)])
    disp('- - - - - -')
end


%% PRINT all subfigs
if exist('print_flag','var')&&(print_flag==1)
    subfigs = fieldnames(figs);
    for i = 1:length(subfigs)
        for j = 1:numel(figs.(subfigs{i}))
            print(figs.(subfigs{i})(j),[savedir,'/aKO_',subfigs{i},'(',num2str(j),').eps'], '-depsc')      
        end
    end
    disp(['saved subfigs to ''',savedir,'aKO_...'''])
end