
%% EMO ANALYSIS 
clear all
cd('/Users/kai/UOL_MPPL/Projects/VibroChair/Vibrotactile_Experiment1/Analysis_Final/AllFiles') % adjust path here 

% get data
import_EMO_exp1
import_EMO_exp2
cohort2.SubjID = cohort2.SubjID + 22; % make compatible for merging with data from exp1
EMOall = vertcat(cohort1, cohort2);
NumP = 41; 

%% sort data
EMOm = grpstats(EMOall, ["SubjID","CondIdx"], 'mean'); % average across excerpts 
vars = {'Pleasantness', 'Arousal', 'Groove', 'LivePerformance', 'BeingPart', 'Liking'}; 
varnames = {'Valence', 'Arousal', 'Groove', 'LiveFeeling', 'BeingPart', 'Liking'}; 

% generate indices
for k = 1:NumP 
    ind.part(:,k) = EMOm.SubjID == k;
end
for k = 1:4
    ind.cond(:,k) = EMOm.CondIdx == k;
end

%% stats modeling

% prepare data
vars = {'Pleasantness', 'Arousal', 'Groove', 'LivePerformance', 'BeingPart', 'Liking'}; 

for nV = 1:6
    ds = table(EMOall.(vars{nV}), categorical(EMOall.SongIdx), categorical(EMOall.CondIdx), categorical(EMOall.SubjID));
    ds.Properties.VariableNames = {'dv', 'song', 'vib', 'part'};
    LME = fitlme(ds, 'dv ~ 1 + vib + (1|part) + (1|song)', 'DummyVarCoding', 'reference'); 
    res.LME(nV).x = LME; 
end

% get coefficients and CIs 
for nV = 1:6
    betas(:,nV) = res.LME(nV).x.Coefficients.Estimate(2:4);
    cilow(:,nV) = res.LME(nV).x.Coefficients.Lower(2:4);
    ciupp(:,nV) = res.LME(nV).x.Coefficients.Upper(2:4);
end

for nV = 1:6
    res.LME(nV).x.Rsquared
end

%% plot emo data coefficients estimates
cmm = colormap('lines'); 
cm = cmm;
cm(1,:) = cmm(1,:).^3; 
cm(2,:) = cmm(1,:).^1; 
cm(3,:) = cmm(1,:).^.25; 
cm(4,:) = [.5 .5 .5]; 

figure; 


% radar/glyph plot!
for nC = 1:4
    for nVar = 1:6
        gplot(nC,nVar) = nanmean(EMOm.(strcat('mean_', vars{nVar}))(ind.cond(:, nC)));
    end
end

subplot(1,2,1);
spider_plot(gplot, 'linewidth', 4, 'color', cm, 'AxesLabels', varnames, 'AxesLabelsEdge', 'none', 'LabelFontSize', 20,...
    'AxesLimits', [2.5*ones(1,6); 5.5*ones(1,6)], 'axesshadedlimits',  {[2.5*ones(1,6); 5.5*ones(1,6)]});
text(-1.4, 1.3, 'A', 'fontsize', 24, 'fontweight', 'bold')

% barplot of betas 
subplot(1,2,2); 
hold on
xpos = -100;
for nC = 1:4
    plot([xpos, xpos], [-1, -4], '-', 'linewidth', 6, 'color', cm(nC, :));
end

for nC = 1:3
    for nVar = 1:6
        m = betas(nC,nVar);
        x = EMOm.(strcat('mean_', vars{nVar}))(ind.cond(:, nC));
        res.mean_g1(nC,nVar) = mean(x(1:22));
        res.mean_g2(nC,nVar) = mean(x(23:end));
        CI = [cilow(nC,nVar), ciupp(nC,nVar)]; 
        xpos = .18*nC+nVar - .2; 
        plot([xpos, xpos], [0, m], '-', 'linewidth', 12, 'color', cm(nC+1, :));
        plot([xpos, xpos], [CI(1), CI(2)], '-', 'linewidth', 3, 'color', [1 1 1]); hold on
        plot([xpos, xpos], [CI(1), CI(2)], '-', 'linewidth', 2, 'color', [cm(nC+1,:)]); hold on
    end
end

ylim([-2.2, 2.2])
xlim([0.5, 6.5])
legend('Headph', 'Mono', 'Multi', 'Incongr')
set(gca, 'XTick', [1:6])
set(gca, 'XTickLabel', varnames)
ff = gca; 
ff.XTickLabelRotation = 30;
ylabel('\beta')
text(-.5, 2.2, 'B', 'fontsize', 24, 'fontweight', 'bold')
subplot(1,2,2)
grid off

%% group/cohort differences? 
corrplot2(squish(res.mean_g1(:,:)), squish(res.mean_g2(:,:)))
[r, p, ci] = corrcoef(squish(res.mean_g1(:,:)), squish(res.mean_g2(:,:)))

%% compute PCA, explore dimensionality
emo = table2array(EMOm);
E = emo(:,5:10);

qf_bartlett(E) % i.e. covariance matrix different than identity matrix 

% look at correlation matrix 
[c, p] = corrcoef(E)
corrr = tril(corrcoef(E), -1); 
cvals = corrr(corrr~=0);

% do PCA 
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca((E)); 
figure; 
subplot(1,2,1)
plot(EXPLAINED)
cumsum(EXPLAINED)

% print coeffs 
COEFF

% varimax rotation 
[fact, rot] = rotatefactors(COEFF(:,1:2), 'Method', 'varimax');
%new_SCORE = SCORE(:,1:2)*rot; % same as new_SCORE = Z*COEFF*rot = Z*fact
%round(fact, 2)
COEFF(:,1:2)*rot
newSc = SCORE(:,1:2)*rot

% print rotated factors 
round(COEFF(:,1:2)*rot, 2)

% reverse sign for ease of interpretability 
for nVib = 1:4
    newSc1(:, nVib) = newSc(ind.cond(:,nVib),1); 
    newSc2(:, nVib) = -newSc(ind.cond(:,nVib),2); % reverse signs here! 
end


%% plot factors 
linsty = {'o', 's', '>', 'd'}; 

figure; 
subplot(1,2,1)
shift_text = [0 .08; 0 0; 0 0; 0 0; 0 0; 0 0];

fac = COEFF(:,1:2)*rot; 
for nFac = 1:6
    plot([0 fac(nFac, 1)], -[0, fac(nFac,2)], 'linewidth', 1, 'color', cm(1,:)); hold on
    plot(fac(nFac, 1), -fac(nFac,2), 'o', 'linewidth', 2, 'color', cm(1,:))
    text(fac(nFac, 1)+.05, -fac(nFac,2)+shift_text(nFac,2), varnames{nFac}, 'fontsize', 16)
end
xlim([-.5 1.2])
ylim([-.5 1.2])
yline(0)
xline(0)
set(gca, 'Box', 'off')
xlabel('Factor 1: Engagement')
ylabel('Factor 2: Preference')
text(-.9, 1.3, 'A', 'fontsize', 24, 'fontweight', 'bold')


subplot(1,2,2)
hold on
m1 = mean(newSc1); % get the data straight 
m2 = mean(newSc2);
ci1 = boot_CI(newSc1);
ci2 = boot_CI(newSc2); 

% plot for legend 
for nVib = 1:4
    plot(m1(nVib), m2(nVib), 'LineStyle', 'none', 'Marker', linsty{nVib}, ...
       'MarkerFaceColor',cm(nVib,:),'MarkerEdgeColor', [.1 .1 .1], ...
       'MarkerSize', 14); hold on
end

for nVib = 1:4
    plot([m1(nVib) m1(nVib)], [ci2(1,nVib), ci2(2,nVib)], '-', 'linewidth', 3, 'color', cm(nVib, :))
    plot([ci1(1,nVib), ci1(2,nVib)], [m2(nVib) m2(nVib)], '-', 'linewidth', 3, 'color', cm(nVib, :))
end
for nVib = 1:4
    plot(m1(nVib), m2(nVib), ...
   'LineStyle', 'none', 'Marker', linsty{nVib}, ...
       'MarkerFaceColor',cm(nVib,:),'MarkerEdgeColor', [.1 .1 .1], ...
       'MarkerSize', 14); hold on
end

xlabel('Factor 1: Engagement')
ylabel('Factor 2: Preference')
legend('Headph', 'Mono', 'Multi', 'Incongr')
axis([-2 2 -2 2])

% interpretation
% dim 1: engagement [1. being part, 2. groove, 3. arousal, 4. live performance]
% dim 2: unpleasantness [1.-pleasantness, 2.-liking, 3.-life performance, 4.-groove]
% [reverse sign for latter]
text(-2.7, 2.2, 'B', 'fontsize', 24, 'fontweight', 'bold')



%% indiv diffs -- gradients in 2d space? 
% import demographics
cd('/Users/kai/UOL_MPPL/Projects/VibroChair/Vibrotactile_Experiment1/Analysis_Final/AllFiles')
import_demographics

res.demo.MSI = mean([demographics.MSI_Train, demographics.MSI_Percept],2); 
res.demo.age = (demographics.Age); 

%% demographics scatterplots 
figure
tittext = {'Headph', 'Mono', 'Multi', 'Incongr'};
for nVib = 1:4
    subplot(2,4, nVib); hold on
    nicescatter(res.demo.MSI(:), newSc1(:,nVib), 'o', cm(nVib,:))
    [r,p] = corrplotXY(res.demo.MSI(:), newSc1(:,nVib))
    title(tittext{nVib})
    xlim([10 50]); ylim([-4 4])
    res.corr.ind_diffs.fac1.r(nVib) = r; 
    res.corr.ind_diffs.fac1.p(nVib) = p; 

    subplot(2,4, 4+nVib); hold on
    nicescatter(res.demo.MSI(:), newSc2(:,nVib), 'o', cm(nVib,:))
    [r,p] = corrplotXY(res.demo.MSI(:), newSc2(:,nVib))
    xlim([10 50]); ylim([-4 4])
    res.corr.ind_diffs.fac2.r(nVib) = r; 
    res.corr.ind_diffs.fac2.p(nVib) = p; 
end

subplot(2,4,1);
ylabel('Engagement')
xlabel('MSI')

subplot(2,4,5);
ylabel('Preference')
xlabel('MSI')

%% post hoc tests for differences of scores?! 
% mono and multi are not different in 2d
[h,p, ci, stats] = ttest(newSc1(:, 2),newSc1(:, 3))
[h,p, ci, stats] = ttest(newSc2(:, 2),newSc2(:, 3))

% headphones and mono are different in 2d
[h,p, ci, stats] = ttest(newSc1(:, 1),newSc1(:, 2))
[h,p, ci, stats] = ttest(newSc2(:, 1),newSc2(:, 2))

% headphones and multi are different in 2d
[h,p, ci, stats] = ttest(newSc1(:, 1),newSc1(:, 3))
[h,p, ci, stats] = ttest(newSc2(:, 1),newSc2(:, 3))

% headphones and incongruent are different in 2d
[h,p, ci, stats] = ttest(newSc1(:, 1),newSc1(:, 4))
[h,p, ci, stats] = ttest(newSc2(:, 1),newSc2(:, 4))


% incongr and mono are different in 2d
[h,p, ci, stats] = ttest(newSc1(:, 4),newSc1(:, 2))
[h,p, ci, stats] = ttest(newSc2(:, 4),newSc2(:, 2))

% incongr and multi are different in 2d
[h,p, ci, stats] = ttest(newSc1(:, 4),newSc1(:, 3))
[h,p, ci, stats] = ttest(newSc2(:, 4),newSc2(:, 3))


