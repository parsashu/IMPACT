function outputs = analyze_amplitude_fidelity(resultsFile, outputDir)
% PACT_PLOT_AMPLITUDE_REVISION_V2
% Reviewer-driven amplitude-fidelity analysis for PACT v9.3.
%
% Reads all_reconstruction_results(_PARTIAL).csv and produces:
%   1) Figure_amplitude_optimal_scaling.png/pdf
%      - NRMSE after per-algorithm optimal zero-intercept global scaling
%      - target-amplitude recovery after the same scaling
%      - full-view and limited-view panels
%   2) Table_amplitude_representative.csv
%      - 20, 0, and -5 dB summary values (phantom-cluster means)
%   3) Table_depth_recovery.csv
%      - inner/middle/outer radial target-amplitude recovery
%   4) Table_common_DAS_scaling.csv
%      - strict DAS-derived common-scale stress test
%
% IMPORTANT INTERPRETATION
% The primary quantitative-amplitude analysis is the per-algorithm optimal
% zero-intercept scale, alpha*=<r,p0>/<r,r>, because nonlinear beamformers
% can have different native output scales/units. The DAS-derived common
% scale is retained as a strict no-recalibration stress test and should be
% reported as supplementary evidence, not as the sole amplitude metric.
%
% Example:
%   outputs = analyze_amplitude_fidelity( ...
%       'all_reconstruction_results_PARTIAL.csv','AMPLITUDE_REVISION');

if nargin < 1 || strlength(string(resultsFile)) == 0
    resultsFile = 'all_reconstruction_results_PARTIAL.csv';
end
if nargin < 2 || strlength(string(outputDir)) == 0
    outputDir = 'AMPLITUDE_REVISION';
end
resultsFile = char(resultsFile);
outputDir = char(outputDir);
if ~isfile(resultsFile), error('Results file not found: %s',resultsFile); end
if ~isfolder(outputDir), mkdir(outputDir); end

T = readtable(resultsFile,'VariableNamingRule','preserve');
T.Phantom = string(T.Phantom);
T.SensorGeometry = string(T.SensorGeometry);
T.Algorithm = string(T.Algorithm);

required = ["Phantom","SensorGeometry","Algorithm","RequestedSNR_dB", ...
    "NoiseRealization","NRMSE_OptimalScale","TargetRecovery_OptimalScale", ...
    "NRMSE_DASCommonScale","TargetRecovery_DASCommonScale", ...
    "TargetRecovery_Inner","TargetRecovery_Middle","TargetRecovery_Outer"];
missing = setdiff(required,string(T.Properties.VariableNames));
if ~isempty(missing)
    error('Missing required columns: %s',strjoin(missing,', '));
end

algorithms = ["DAS","CF-DAS","SCF-DAS","DMAS","DASDSF","NC-CF","NC-VDAS"];
T = T(ismember(T.Algorithm,algorithms),:);

% First average repeated noise realizations within each phantom/condition.
cluster = phantomClusterMeans(T);
summary = summarizeAcrossPhantoms(cluster,4000,661177);

% Representative table requested by Reviewer 2.
representativeSNR = [20,0,-5];
rep = summary(ismember(summary.RequestedSNR_dB,representativeSNR),:);
rep = sortrows(rep,{'SensorGeometry','RequestedSNR_dB','Algorithm'});
writetable(rep,fullfile(outputDir,'Table_amplitude_representative.csv'));

% Depth/radial-position recovery table. Inner corresponds to deeper/central
% targets; outer corresponds to more superficial/peripheral targets.
depthVars = {'Phantom','SensorGeometry','RequestedSNR_dB','Algorithm', ...
    'TargetRecovery_Inner','TargetRecovery_Middle','TargetRecovery_Outer'};
depth = cluster(:,depthVars);
writetable(depth,fullfile(outputDir,'Table_depth_recovery.csv'));

% Strict common-DAS scaling table (supplementary/stress-test use).
commonVars = {'Phantom','SensorGeometry','RequestedSNR_dB','Algorithm', ...
    'NRMSE_DASCommonScale','TargetRecovery_DASCommonScale'};
common = cluster(:,commonVars);
writetable(common,fullfile(outputDir,'Table_common_DAS_scaling.csv'));

% Publication figure: 4 panels. Values are means over phantom-cluster means;
% error bars are approximate 95% CIs across phantoms.
fig = figure('Color','w','Position',[60 40 1700 1050]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
plotPanel(nexttile(tl,1),summary,"CIRCULAR_FULL",algorithms, ...
    'NRMSE','(a) Full-view: amplitude NRMSE');
plotPanel(nexttile(tl,2),summary,"LIMITED_VIEW_ARC",algorithms, ...
    'NRMSE','(b) Limited-view: amplitude NRMSE');
plotPanel(nexttile(tl,3),summary,"CIRCULAR_FULL",algorithms, ...
    'RECOVERY','(c) Full-view: target-amplitude recovery');
plotPanel(nexttile(tl,4),summary,"LIMITED_VIEW_ARC",algorithms, ...
    'RECOVERY','(d) Limited-view: target-amplitude recovery');
title(tl,['Amplitude-sensitive evaluation before independent image normalization' newline ...
    'Optimal zero-intercept global scaling; phantom-cluster mean and 95% bootstrap CI'], ...
    'FontWeight','bold');

pngFile = fullfile(outputDir,'Figure_amplitude_optimal_scaling.png');
pdfFile = fullfile(outputDir,'Figure_amplitude_optimal_scaling.pdf');
exportgraphics(fig,pngFile,'Resolution',300);
exportgraphics(fig,pdfFile,'ContentType','vector');
close(fig);

outputs = struct;
outputs.summary = summary;
outputs.representative = rep;
outputs.depth = depth;
outputs.common = common;
outputs.figurePNG = string(pngFile);
outputs.figurePDF = string(pdfFile);
save(fullfile(outputDir,'amplitude_revision_outputs.mat'),'outputs','-v7.3');

fprintf('\nAmplitude revision outputs saved to:\n  %s\n',outputDir);
end

function C = phantomClusterMeans(T)
keys = {'Phantom','SensorGeometry','RequestedSNR_dB','Algorithm'};
[G,p,g,s,a] = findgroups(T.Phantom,T.SensorGeometry,T.RequestedSNR_dB,T.Algorithm);
C = table(p,g,s,a,'VariableNames',keys);
metrics = ["NRMSE_OptimalScale","TargetRecovery_OptimalScale", ...
    "NRMSE_DASCommonScale","TargetRecovery_DASCommonScale", ...
    "TargetRecovery_Inner","TargetRecovery_Middle","TargetRecovery_Outer"];
for k = 1:numel(metrics)
    C.(metrics(k)) = splitapply(@meanFinite,T.(metrics(k)),G);
end
end

function S = summarizeAcrossPhantoms(C,bootstrapCount,statisticsSeed)
[G,g,s,a] = findgroups(C.SensorGeometry,C.RequestedSNR_dB,C.Algorithm);
S = table(g,s,a,'VariableNames',{'SensorGeometry','RequestedSNR_dB','Algorithm'});
metrics = ["NRMSE_OptimalScale","TargetRecovery_OptimalScale", ...
    "NRMSE_DASCommonScale","TargetRecovery_DASCommonScale"];

% Keep inference at the phantom level. Noise realizations have already been
% averaged in phantomClusterMeans. Bootstrap the phantom-level observations
% to match the manuscript's inferential framework.
rng(statisticsSeed,'twister');
for k = 1:numel(metrics)
    x = C.(metrics(k));
    S.(metrics(k)+"_Mean") = splitapply(@meanFinite,x,G);
    S.(metrics(k)+"_SD") = splitapply(@stdFinite,x,G);
    S.(metrics(k)+"_N") = splitapply(@countFinite,x,G);
    lo = nan(height(S),1);
    hi = nan(height(S),1);
    for j = 1:height(S)
        vals = x(G==j);
        vals = vals(isfinite(vals));
        n = numel(vals);
        if n == 0
            continue;
        elseif n == 1
            lo(j) = vals(1);
            hi(j) = vals(1);
            continue;
        end
        idx = randi(n,n,bootstrapCount);
        bootMeans = mean(vals(idx),1);
        lo(j) = percentileLinearLocal(bootMeans,2.5);
        hi(j) = percentileLinearLocal(bootMeans,97.5);
    end
    S.(metrics(k)+"_CI95Low") = lo;
    S.(metrics(k)+"_CI95High") = hi;
end
end

function plotPanel(ax,S,geometry,algorithms,metricType,panelTitle)
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
% Display SNR from high-quality/noiseless to severe noise.
snrOrder = [Inf,30,20,10,5,0,-5];
x = 1:numel(snrOrder);
for k = 1:numel(algorithms)
    alg = algorithms(k);
    mu = nan(size(x)); lo = nan(size(x)); hi = nan(size(x));
    for q = 1:numel(snrOrder)
        idx = S.SensorGeometry==geometry & S.Algorithm==alg & ...
            snrEqualVector(S.RequestedSNR_dB,snrOrder(q));
        if ~any(idx), continue; end
        if metricType=="NRMSE"
            mu(q) = S.NRMSE_OptimalScale_Mean(idx);
            lo(q) = S.NRMSE_OptimalScale_CI95Low(idx);
            hi(q) = S.NRMSE_OptimalScale_CI95High(idx);
        else
            mu(q) = S.TargetRecovery_OptimalScale_Mean(idx);
            lo(q) = S.TargetRecovery_OptimalScale_CI95Low(idx);
            hi(q) = S.TargetRecovery_OptimalScale_CI95High(idx);
        end
    end
    e = errorbar(ax,x,mu,mu-lo,hi-mu,'-o','LineWidth',1.25, ...
        'MarkerSize',4,'DisplayName',char(alg));
    if alg=="DAS", e.LineWidth = 2.0; end
end
xticks(ax,x); xticklabels(ax,{'Inf','30','20','10','5','0','-5'});
xlabel(ax,'Input SNR (dB)');
if metricType=="NRMSE"
    ylabel(ax,'NRMSE after optimal global scaling');
else
    ylabel(ax,'Relative target-amplitude recovery');
    yline(ax,1,'--','Ideal recovery = 1','HandleVisibility','off');
    ylim(ax,[0 1.05]);
end
title(ax,panelTitle,'FontWeight','bold');
if contains(panelTitle,'(a)')
    legend(ax,'Location','eastoutside');
end
end

function tf = snrEqualVector(values,target)
if isinf(target)
    tf = isinf(values) & sign(values)==sign(target);
else
    tf = abs(values-target)<1e-9;
end
end

function y = meanFinite(x)
x = x(isfinite(x));
if isempty(x), y = NaN; else, y = mean(x); end
end

function y = stdFinite(x)
x = x(isfinite(x));
if numel(x)<2, y = NaN; else, y = std(x,0); end
end

function y = countFinite(x)
y = nnz(isfinite(x));
end

function q = percentileLinearLocal(x,p)
x = sort(double(x(:)));
if isempty(x)
    q = NaN;
    return;
elseif numel(x)==1
    q = x(1);
    return;
end
pos = 1 + (numel(x)-1)*(p/100);
lo = floor(pos);
hi = ceil(pos);
if lo==hi
    q = x(lo);
else
    w = pos-lo;
    q = (1-w)*x(lo) + w*x(hi);
end
end
