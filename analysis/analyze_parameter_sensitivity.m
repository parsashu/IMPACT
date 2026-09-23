function outputs = analyze_parameter_sensitivity( ...
    sensitivityCsv, mainResultsCsv, outputDir)
% analyze_parameter_sensitivity
% Reviewer-driven parameter-sensitivity post-processing.
%
% The analysis averages repeated noise realizations first within each
% phantom/geometry/SNR/parameter setting, then summarizes the resulting
% phantom-level condition means. Matching DAS rows from the main benchmark
% are used only to compute metric differences versus DAS.
%
% Required inputs:
%   parameter_sensitivity_PARTIAL.csv
%   all_reconstruction_results_PARTIAL.csv
%
% The "_PARTIAL" sensitivity file is complete for the configured sweep:
% 1593 rows + header.

if nargin < 1 || isempty(sensitivityCsv)
    sensitivityCsv = 'parameter_sensitivity_PARTIAL.csv';
end
if nargin < 2 || isempty(mainResultsCsv)
    mainResultsCsv = 'all_reconstruction_results_PARTIAL.csv';
end
if nargin < 3 || isempty(outputDir)
    outputDir = 'PARAMETER_SENSITIVITY_REVISION_MATLAB';
end
if ~isfolder(outputDir), mkdir(outputDir); end

S = readtable(sensitivityCsv,'TextType','string');
M = readtable(mainResultsCsv,'TextType','string');

requiredS = ["Phantom","SensorGeometry","RequestedSNR_dB","NoiseRealization", ...
    "ParameterFamily","Parameter1Name","Parameter1Value","Parameter2Name", ...
    "Parameter2Value","Algorithm","PSNR_ROI_dB","SSIM","Dice", ...
    "ArtifactToTarget_dB","Runtime_Sec"];
assert(all(ismember(requiredS,string(S.Properties.VariableNames))), ...
    'Sensitivity CSV is missing required columns.');

% Matching DAS reference cases.
phantoms = unique(S.Phantom,'stable');
snrs = unique(S.RequestedSNR_dB,'stable');
reps = unique(S.NoiseRealization,'stable');
idxD = M.Algorithm=="DAS" & ismember(M.Phantom,phantoms) & ...
    ismember(M.RequestedSNR_dB,snrs) & ismember(M.NoiseRealization,reps);
DAS = M(idxD,["Phantom","SensorGeometry","RequestedSNR_dB", ...
    "NoiseRealization","PSNR_ROI_dB","SSIM","Dice","ArtifactToTarget_dB"]);
DAS.Properties.VariableNames(end-3:end) = ...
    {'DAS_PSNR_ROI_dB','DAS_SSIM','DAS_Dice','DAS_ArtifactToTarget_dB'};

T = outerjoin(S,DAS,'Keys', ...
    {'Phantom','SensorGeometry','RequestedSNR_dB','NoiseRealization'}, ...
    'MergeKeys',true,'Type','left');

T.DeltaPSNR_vs_DAS_dB = T.PSNR_ROI_dB - T.DAS_PSNR_ROI_dB;
T.DeltaSSIM_vs_DAS = T.SSIM - T.DAS_SSIM;
T.DeltaDice_vs_DAS = T.Dice - T.DAS_Dice;
T.DeltaArtifact_vs_DAS_dB = ...
    T.ArtifactToTarget_dB - T.DAS_ArtifactToTarget_dB;

% Average realizations first.
G = findgroups(T.Phantom,T.SensorGeometry,T.RequestedSNR_dB, ...
    T.ParameterFamily,T.Parameter1Name,T.Parameter1Value, ...
    T.Parameter2Name,T.Parameter2Value,T.Algorithm);

P = table();
[P.Phantom,P.SensorGeometry,P.RequestedSNR_dB,P.ParameterFamily, ...
 P.Parameter1Name,P.Parameter1Value,P.Parameter2Name,P.Parameter2Value, ...
 P.Algorithm] = splitapply(@firstTuple, ...
    T.Phantom,T.SensorGeometry,T.RequestedSNR_dB,T.ParameterFamily, ...
    T.Parameter1Name,T.Parameter1Value,T.Parameter2Name,T.Parameter2Value, ...
    T.Algorithm,G);

metricNames = ["PSNR_ROI_dB","SSIM","Dice","ArtifactToTarget_dB","Runtime_Sec", ...
    "DeltaPSNR_vs_DAS_dB","DeltaSSIM_vs_DAS","DeltaDice_vs_DAS", ...
    "DeltaArtifact_vs_DAS_dB"];
for k = 1:numel(metricNames)
    P.(metricNames(k)) = splitapply(@(x)mean(x,'omitnan'),T.(metricNames(k)),G);
end
writetable(P,fullfile(outputDir,'parameter_sensitivity_phantom_condition_means.csv'));

% Setting-level summary.
G2 = findgroups(P.ParameterFamily,P.Parameter1Name,P.Parameter1Value, ...
    P.Parameter2Name,P.Parameter2Value,P.Algorithm);
Q = table();
[Q.ParameterFamily,Q.Parameter1Name,Q.Parameter1Value, ...
 Q.Parameter2Name,Q.Parameter2Value,Q.Algorithm] = ...
    splitapply(@firstTuple,P.ParameterFamily,P.Parameter1Name, ...
    P.Parameter1Value,P.Parameter2Name,P.Parameter2Value,P.Algorithm,G2);
Q.NumClusters = splitapply(@numel,P.Phantom,G2);
for k = 1:numel(metricNames)
    Q.("Mean_"+metricNames(k)) = ...
        splitapply(@(x)mean(x,'omitnan'),P.(metricNames(k)),G2);
end
Q.PrimarySetting = false(height(Q),1);
Q.PrimarySetting(Q.ParameterFamily=="SCF_EXPONENT" & Q.Parameter1Value==1.5) = true;
Q.PrimarySetting(Q.ParameterFamily=="DASDSF_EXPONENT" & Q.Parameter1Value==1.0) = true;
Q.PrimarySetting(Q.ParameterFamily=="NCVDAS_PARAMETERS" & ...
    Q.Parameter1Value==1.0 & Q.Parameter2Value==1.0) = true;
Q.PrimarySetting(Q.ParameterFamily=="MV_PARAMETERS" & ...
    Q.Parameter1Value==16 & abs(Q.Parameter2Value-0.01)<eps) = true;
writetable(Q,fullfile(outputDir,'parameter_sensitivity_setting_summary.csv'));

% Geometry/SNR summary for delta metrics.
G3 = findgroups(P.ParameterFamily,P.Parameter1Value,P.Parameter2Value, ...
    P.SensorGeometry,P.RequestedSNR_dB);
R = table();
[R.ParameterFamily,R.Parameter1Value,R.Parameter2Value, ...
 R.SensorGeometry,R.RequestedSNR_dB] = ...
    splitapply(@firstTuple,P.ParameterFamily,P.Parameter1Value, ...
    P.Parameter2Value,P.SensorGeometry,P.RequestedSNR_dB,G3);
R.MeanDeltaPSNR_vs_DAS_dB = splitapply(@(x)mean(x,'omitnan'), ...
    P.DeltaPSNR_vs_DAS_dB,G3);
R.MeanDeltaSSIM_vs_DAS = splitapply(@(x)mean(x,'omitnan'), ...
    P.DeltaSSIM_vs_DAS,G3);
R.MeanDeltaDice_vs_DAS = splitapply(@(x)mean(x,'omitnan'), ...
    P.DeltaDice_vs_DAS,G3);
writetable(R,fullfile(outputDir,'parameter_sensitivity_vs_DAS_by_geometry_SNR.csv'));

% -------------------------------------------------------------------------
% Figure 1: SCF
fig = figure('Color','w','Position',[50 50 900 600]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plotOneDimSensitivity(ax,P,R,"SCF_EXPONENT",1.5);
xlabel(ax,'SCF exponent p'); ylabel(ax,'Mean \DeltaSSIM relative to DAS');
title(ax,'SCF parameter sensitivity');
legend(ax,'Location','best','Box','off');
exportgraphics(fig,fullfile(outputDir,'Fig_parameter_sensitivity_SCF.png'),'Resolution',300);
exportgraphics(fig,fullfile(outputDir,'Fig_parameter_sensitivity_SCF.pdf'),'ContentType','vector');
close(fig);

% Figure 2: DASDSF
fig = figure('Color','w','Position',[50 50 900 600]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plotOneDimSensitivity(ax,P,R,"DASDSF_EXPONENT",1.0);
xlabel(ax,'DASDSF exponent q'); ylabel(ax,'Mean \DeltaSSIM relative to DAS');
title(ax,'DASDSF parameter sensitivity');
legend(ax,'Location','best','Box','off');
exportgraphics(fig,fullfile(outputDir,'Fig_parameter_sensitivity_DASDSF.png'),'Resolution',300);
exportgraphics(fig,fullfile(outputDir,'Fig_parameter_sensitivity_DASDSF.pdf'),'ContentType','vector');
close(fig);

% Figure 3: NC-VDAS heatmap
makeHeatmap(Q,"NCVDAS_PARAMETERS", ...
    'Noise-penalty coefficient \beta','NC-VDAS exponent q', ...
    'NC-VDAS: mean \DeltaSSIM relative to DAS', ...
    [1 1],fullfile(outputDir,'Fig_parameter_sensitivity_NCVDAS.png'), ...
    fullfile(outputDir,'Fig_parameter_sensitivity_NCVDAS.pdf'));

% Figure 4: MV heatmap
makeHeatmap(Q,"MV_PARAMETERS", ...
    'Diagonal loading \delta','Subarray length L', ...
    'SM-MV: mean \DeltaSSIM relative to DAS (limited view)', ...
    [16 0.01],fullfile(outputDir,'Fig_parameter_sensitivity_MV.png'), ...
    fullfile(outputDir,'Fig_parameter_sensitivity_MV.pdf'));

outputs = struct('phantomConditionMeans',P,'settingSummary',Q, ...
    'geometrySNRSummary',R,'outputDir',string(outputDir));
end

function varargout = firstTuple(varargin)
varargout = cell(size(varargin));
for k = 1:nargin
    varargout{k} = varargin{k}(1);
end
end

function plotOneDimSensitivity(ax,P,R,family,primaryValue)
% Overall across all phantom-geometry-SNR clusters
idxP = P.ParameterFamily==family;
values = unique(P.Parameter1Value(idxP));
values = sort(values);
overall = nan(size(values));
fullm5 = nan(size(values));
limm5 = nan(size(values));
for k = 1:numel(values)
    overall(k) = mean(P.DeltaSSIM_vs_DAS(idxP & P.Parameter1Value==values(k)),'omitnan');
    idx = R.ParameterFamily==family & R.Parameter1Value==values(k) & ...
        R.SensorGeometry=="CIRCULAR_FULL" & R.RequestedSNR_dB==-5;
    fullm5(k) = R.MeanDeltaSSIM_vs_DAS(idx);
    idx = R.ParameterFamily==family & R.Parameter1Value==values(k) & ...
        R.SensorGeometry=="LIMITED_VIEW_ARC" & R.RequestedSNR_dB==-5;
    limm5(k) = R.MeanDeltaSSIM_vs_DAS(idx);
end
plot(ax,values,overall,'-o','LineWidth',1.5,'DisplayName','All conditions');
plot(ax,values,fullm5,'-o','LineWidth',1.5,'DisplayName','Full view, -5 dB');
plot(ax,values,limm5,'-o','LineWidth',1.5,'DisplayName','Limited view, -5 dB');
yline(ax,0,'--','HandleVisibility','off');
xline(ax,primaryValue,':','Primary setting','HandleVisibility','off');
end

function makeHeatmap(Q,family,xlab,ylab,ttl,primary,pngFile,pdfFile)
S = Q(Q.ParameterFamily==family,:);
xv = sort(unique(S.Parameter2Value));
yv = sort(unique(S.Parameter1Value));
Z = nan(numel(yv),numel(xv));
for i=1:numel(yv)
    for j=1:numel(xv)
        idx = S.Parameter1Value==yv(i) & S.Parameter2Value==xv(j);
        Z(i,j) = S.Mean_DeltaSSIM_vs_DAS(idx);
    end
end
fig = figure('Color','w','Position',[50 50 800 600]);
h = heatmap(string(xv),string(yv),Z);
h.XLabel = xlab; h.YLabel = ylab; h.Title = ttl;
h.CellLabelFormat = '%+.4f';
exportgraphics(fig,pngFile,'Resolution',300);
exportgraphics(fig,pdfFile,'ContentType','vector');
close(fig);
end
