extracted_2state_stats_file = 'extracted_2state_stats.mat';

if ~exist(extracted_2state_stats_file)
    % Load 2-state params
    sourcedir = './results/';
    ncells = 28;
    D_bin = nan(ncells,2); % diffusivity
    D_tag = nan(ncells,2);
    P_bin = nan(ncells,2); % occupancy
    P_tag = nan(ncells,2);
    T_bin = nan(ncells,2); % dwell time
    T_tag = nan(ncells,2);
    A_bin = nan(ncells,2,2); % transition matrix
    A_tag = nan(ncells,2,2);

    for cellid=1:ncells
        fprintf('Extracting 2-state model statistics for cell %i\n',cellid);
        f_bin=dir([sourcedir sprintf('bin_%02d',cellid) '_HMManalysis_hidden2.mat']);
        f_tag=dir([sourcedir sprintf('tag_%02d',cellid) '_HMManalysis_hidden2.mat']);
        files_bin = {f_bin.name}';
        files_tag = {f_tag.name}';

        result_tag = load([sourcedir sprintf('tag_%02d_HMManalysis_hidden2.mat',cellid)]);
        result_bin = load([sourcedir sprintf('bin_%02d_HMManalysis_hidden2.mat',cellid)]);

        D_tag(cellid,:) = result_tag.WbestN{2}.est.DdtMean/0.02;
        D_bin(cellid,:) = result_bin.WbestN{2}.est.DdtMean/0.02;

        P_tag(cellid,:) = result_tag.WbestN{2}.est.Ptot;
        P_bin(cellid,:) = result_bin.WbestN{2}.est.Ptot;
        
        T_tag(cellid,:) = result_tag.WbestN{2}.est.dwellMean;
        T_bin(cellid,:) = result_bin.WbestN{2}.est.dwellMean;

        A_tag(cellid,:,:) = result_tag.WbestN{2}.est.Amean;
        A_bin(cellid,:,:) = result_bin.WbestN{2}.est.Amean;
    end
    save(extracted_2state_stats_file,'D_tag','D_bin','P_tag','P_bin','T_tag','T_bin','A_tag','A_bin');
else
    load(extracted_2state_stats_file)
end

figuredir = 'figures/2state/';

mkdir(figuredir)

%% Diffusivity 
figure('position',[ 383 1023 497 240])
subplot(1,2,1); hold on;
bar(1,mean(D_tag(:,1)),'facecolor',[0 .447 .741],'edgecolor','none');
bar(2,mean(D_bin(:,1)),'facecolor',[1 0 0],'edgecolor','none');
errorbar(1,mean(D_tag(:,1)),std(D_tag(:,1)),'linestyle','none','color','k','linewidth',2);
errorbar(2,mean(D_bin(:,1)),std(D_bin(:,1)),'linestyle','none','color','k','linewidth',2);
ylabel('D (\mum^2/s)');
set(gca,'xtick',1.5,'xticklabel','Slow state')
set(gca,'fontsize',18,'fontname','arial')
set(gca,'ytick',[0 0.01 0.02],'ylim',[0 0.022]);
set(gca,'linewidth',2,'tickdir','out')
subplot(1,2,2); hold on;
bar(1,mean(D_tag(:,2)),'facecolor',[0 .447 .741],'edgecolor','none');
bar(2,mean(D_bin(:,2)),'facecolor',[1 0 0],'edgecolor','none');
errorbar(1,mean(D_tag(:,2)),std(D_tag(:,2)),'linestyle','none','color','k','linewidth',2);
errorbar(2,mean(D_bin(:,2)),std(D_bin(:,2)),'linestyle','none','color','k','linewidth',2);
ylabel('D (\mum^2/s)');
set(gca,'xtick',1.5,'xticklabel','Fast state')
set(gca,'fontsize',18,'fontname','arial')
set(gca,'ytick',[0 0.05 0.10 0.15],'yticklabel',{'0','0.05','0.10','0.15'},'ylim',[0 0.17]);
set(gca,'linewidth',2,'tickdir','out')
filebase = [figuredir 'diffusivity'];
savefig(filebase); % Save the .fig file
print(gcf,'-dtiff',[filebase '.tiff'],'-r300'); % Save a 300 dpi tiff file

%% Occupancy
figure('position',[ 383 1023 497 240])
subplot(1,2,1); hold on;
bar(1,mean(P_tag(:,1)*100),'facecolor',[0 .447 .741],'edgecolor','none');
bar(2,mean(P_bin(:,1)*100),'facecolor',[1 0 0],'edgecolor','none');
errorbar(1,mean(P_tag(:,1)*100),std(P_tag(:,1)*100),'linestyle','none','color','k','linewidth',2);
errorbar(2,mean(P_bin(:,1)*100),std(P_bin(:,1)*100),'linestyle','none','color','k','linewidth',2);
ylabel('Occupancy (%)');
set(gca,'xtick',1.5,'xticklabel','Slow state')
set(gca,'fontsize',18,'fontname','arial')
%set(gca,'ytick',[0 0.01 0.02],'ylim',[0 0.022]);
set(gca,'linewidth',2,'tickdir','out')
subplot(1,2,2); hold on;
bar(1,mean(P_tag(:,2)*100),'facecolor',[0 .447 .741],'edgecolor','none');
bar(2,mean(P_bin(:,2)*100),'facecolor',[1 0 0],'edgecolor','none');
errorbar(1,mean(P_tag(:,2)*100),std(P_tag(:,2)*100),'linestyle','none','color','k','linewidth',2);
errorbar(2,mean(P_bin(:,2)*100),std(P_bin(:,2)*100),'linestyle','none','color','k','linewidth',2);
ylabel('Occupancy (%)');
set(gca,'xtick',1.5,'xticklabel','Fast state')
set(gca,'fontsize',18,'fontname','arial')
%set(gca,'ytick',[0 0.05 0.10 0.15],'yticklabel',{'0','0.05','0.10','0.15'},'ylim',[0 0.17]);
set(gca,'linewidth',2,'tickdir','out')

filebase = [figuredir 'occupancy'];
savefig(filebase); % Save the .fig file
print(gcf,'-dtiff',[filebase '.tiff'],'-r300'); % Save a 300 dpi tiff file

%% Dwell time 
ms_per_frame = 20;
figure('position',[ 383 1023 497 240])
subplot(1,2,1); hold on;
bar(1,mean(T_tag(:,1)*ms_per_frame),'facecolor',[0 .447 .741],'edgecolor','none');
bar(2,mean(T_bin(:,1)*ms_per_frame),'facecolor',[1 0 0],'edgecolor','none');
errorbar(1,mean(T_tag(:,1)*ms_per_frame),std(T_tag(:,1)*ms_per_frame),'linestyle','none','color','k','linewidth',2);
errorbar(2,mean(T_bin(:,1)*ms_per_frame),std(T_bin(:,1)*ms_per_frame),'linestyle','none','color','k','linewidth',2);
ylabel('Dwell time (ms)');
set(gca,'xtick',1.5,'xticklabel','Slow state')
set(gca,'fontsize',18,'fontname','arial')
%set(gca,'ytick',[0 0.01 0.02],'ylim',[0 0.022]);
set(gca,'linewidth',2,'tickdir','out')
subplot(1,2,2); hold on;
bar(1,mean(T_tag(:,2)*ms_per_frame),'facecolor',[0 .447 .741],'edgecolor','none');
bar(2,mean(T_bin(:,2)*ms_per_frame),'facecolor',[1 0 0],'edgecolor','none');
errorbar(1,mean(T_tag(:,2)*ms_per_frame),std(T_tag(:,2)*ms_per_frame),'linestyle','none','color','k','linewidth',2);
errorbar(2,mean(T_bin(:,2)*ms_per_frame),std(T_bin(:,2)*ms_per_frame),'linestyle','none','color','k','linewidth',2);
ylabel('Dwell time (ms)');
set(gca,'xtick',1.5,'xticklabel','Fast state')
set(gca,'fontsize',18,'fontname','arial')
%set(gca,'ytick',[0 0.05 0.10 0.15],'yticklabel',{'0','0.05','0.10','0.15'},'ylim',[0 0.17]);
set(gca,'linewidth',2,'tickdir','out')

filebase = [figuredir 'dwelltime'];
savefig(filebase); % Save the .fig file
print(gcf,'-dtiff',[filebase '.tiff'],'-r300'); % Save a 300 dpi tiff file
