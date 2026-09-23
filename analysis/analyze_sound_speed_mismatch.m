function outputs = analyze_sound_speed_mismatch(inputFile, outputDir)
% PACT_PLOT_DELAY_MISMATCH_REVISION
% Reviewer-driven post-processing for acoustic sound-speed / delay-model
% mismatch. No forward simulation or reconstruction is repeated.
%
% Input: delay_mismatch_results.csv produced by PACT_SoundSpeedMismatch
% Outputs:
%   delay_mismatch_summary.csv
%   delay_mismatch_delta_vs_baseline.csv
%   delay_mismatch_delta_summary.csv
%   delay_mismatch_rank_stability.csv
%   delay_mismatch_perturbation_magnitude.csv
%   figure_delay_mismatch_direct_delta_psnr.pdf/png
%   figure_delay_mismatch_adaptive_limited_psnr.pdf/png
%
% Statistical unit: phantom. Noise realizations are averaged within each
% phantom/condition before across-phantom summaries and bootstrap CIs.

if nargin < 1 || strlength(string(inputFile))==0
    inputFile = "delay_mismatch_results.csv";
end
if nargin < 2 || strlength(string(outputDir))==0
    outputDir = "DELAY_MISMATCH_REVISION";
end
inputFile = string(inputFile);
outputDir = string(outputDir);
if ~isfile(inputFile)
    error('Input file not found: %s',inputFile);
end
if ~isfolder(outputDir), mkdir(outputDir); end

T = readtable(inputFile,'VariableNamingRule','preserve');
required = ["Phantom","MismatchCondition","SensorGeometry","Algorithm", ...
    "RequestedSNR_dB","NoiseRealization","PSNR_ROI_dB","SSIM", ...
    "HD95_mm","gCNR","ArtifactToTarget_dB"];
for k=1:numel(required)
    if ~ismember(required(k),string(T.Properties.VariableNames))
        error('Missing required column: %s',required(k));
    end
end
T.Phantom = string(T.Phantom);
T.MismatchCondition = string(T.MismatchCondition);
T.SensorGeometry = string(T.SensorGeometry);
T.Algorithm = string(T.Algorithm);

% Average repeated noise realizations within each phantom first.
[G,phantom,condition,geometry,snr,algorithm] = findgroups( ...
    T.Phantom,T.MismatchCondition,T.SensorGeometry,T.RequestedSNR_dB,T.Algorithm);
P = table(phantom,condition,geometry,snr,algorithm, ...
    'VariableNames',{'Phantom','MismatchCondition','SensorGeometry', ...
    'RequestedSNR_dB','Algorithm'});
P.PSNR_ROI_dB = splitapply(@meanFinite,T.PSNR_ROI_dB,G);
P.SSIM = splitapply(@meanFinite,T.SSIM,G);
P.HD95_mm = splitapply(@meanFinite,T.HD95_mm,G);
P.gCNR = splitapply(@meanFinite,T.gCNR,G);
P.ArtifactToTarget_dB = splitapply(@meanFinite,T.ArtifactToTarget_dB,G);

% Across-phantom summaries and 95% bootstrap CIs.
rng(661177,'twister');
summary = summarizePhantomMeans(P,4000);
writetable(summary,fullfile(outputDir,'delay_mismatch_summary.csv'));

% Paired differences versus known-map straight-ray baseline.
baseline = P(P.MismatchCondition=="BASELINE_HETERO",:);
base = baseline(:,{'Phantom','SensorGeometry','RequestedSNR_dB','Algorithm', ...
    'PSNR_ROI_dB','SSIM','HD95_mm','gCNR','ArtifactToTarget_dB'});
base.Properties.VariableNames(end-4:end) = ["BasePSNR","BaseSSIM","BaseHD95", ...
    "BasegCNR","BaseArtifactToTarget"];
D = innerjoin(P,base,'Keys',{'Phantom','SensorGeometry','RequestedSNR_dB','Algorithm'});
D.DeltaPSNR_dB = D.PSNR_ROI_dB-D.BasePSNR;
D.DeltaSSIM = D.SSIM-D.BaseSSIM;
D.DeltaHD95_mm = D.HD95_mm-D.BaseHD95;
D.DeltagCNR = D.gCNR-D.BasegCNR;
D.DeltaArtifactToTarget_dB = D.ArtifactToTarget_dB-D.BaseArtifactToTarget;
D = D(D.MismatchCondition~="BASELINE_HETERO",:);
writetable(D,fullfile(outputDir,'delay_mismatch_delta_vs_baseline.csv'));

deltaSummary = summarizeDeltaPhantoms(D,4000);
writetable(deltaSummary,fullfile(outputDir,'delay_mismatch_delta_summary.csv'));

rankStability = computeRankStability(P);
writetable(rankStability,fullfile(outputDir,'delay_mismatch_rank_stability.csv'));

magnitude = makePerturbationMagnitudeTable(T);
writetable(magnitude,fullfile(outputDir,'delay_mismatch_perturbation_magnitude.csv'));

fig1 = drawDirectDeltaPSNR(deltaSummary,outputDir);
fig2 = drawAdaptiveLimitedPSNR(summary,outputDir);

outputs = struct;
outputs.summary = summary;
outputs.delta = D;
outputs.deltaSummary = deltaSummary;
outputs.rankStability = rankStability;
outputs.perturbationMagnitude = magnitude;
outputs.directFigure = fig1;
outputs.adaptiveFigure = fig2;
save(fullfile(outputDir,'delay_mismatch_analysis_outputs.mat'),'outputs','-v7.3');
fprintf('\nDelay-mismatch post-processing complete.\nOutput folder: %s\n',outputDir);
end

function S = summarizePhantomMeans(P,B)
[G,condition,geometry,snr,algorithm] = findgroups(P.MismatchCondition, ...
    P.SensorGeometry,P.RequestedSNR_dB,P.Algorithm);
S = table(condition,geometry,snr,algorithm,'VariableNames', ...
    {'MismatchCondition','SensorGeometry','RequestedSNR_dB','Algorithm'});
metrics = ["PSNR_ROI_dB","SSIM","HD95_mm","gCNR","ArtifactToTarget_dB"];
numGroups=height(S);
for m=1:numel(metrics)
    x=P.(metrics(m));
    means=nan(numGroups,1); sd=nan(numGroups,1); n=nan(numGroups,1);
    lo=nan(numGroups,1); hi=nan(numGroups,1);
    for gg=1:numGroups
        v=x(G==gg); v=v(isfinite(v));
        n(gg)=numel(v);
        if ~isempty(v)
            means(gg)=mean(v);
            if numel(v)>1, sd(gg)=std(v,0); end
            ci=bootstrapPair(v,B); lo(gg)=ci(1); hi(gg)=ci(2);
        end
    end
    S.(metrics(m)+"_Mean")=means;
    S.(metrics(m)+"_SD")=sd;
    S.(metrics(m)+"_N")=n;
    S.(metrics(m)+"_CI95Low")=lo;
    S.(metrics(m)+"_CI95High")=hi;
end
S=sortrows(S,{'SensorGeometry','RequestedSNR_dB','MismatchCondition','Algorithm'});
end

function S = summarizeDeltaPhantoms(D,B)
[G,condition,geometry,snr,algorithm] = findgroups(D.MismatchCondition, ...
    D.SensorGeometry,D.RequestedSNR_dB,D.Algorithm);
S=table(condition,geometry,snr,algorithm,'VariableNames', ...
    {'MismatchCondition','SensorGeometry','RequestedSNR_dB','Algorithm'});
metrics=["DeltaPSNR_dB","DeltaSSIM","DeltaHD95_mm","DeltagCNR","DeltaArtifactToTarget_dB"];
numGroups=height(S);
for m=1:numel(metrics)
    x=D.(metrics(m));
    means=nan(numGroups,1); sd=nan(numGroups,1); n=nan(numGroups,1);
    lo=nan(numGroups,1); hi=nan(numGroups,1);
    for gg=1:numGroups
        v=x(G==gg); v=v(isfinite(v));
        n(gg)=numel(v);
        if ~isempty(v)
            means(gg)=mean(v);
            if numel(v)>1, sd(gg)=std(v,0); end
            ci=bootstrapPair(v,B); lo(gg)=ci(1); hi(gg)=ci(2);
        end
    end
    S.(metrics(m)+"_Mean")=means;
    S.(metrics(m)+"_SD")=sd;
    S.(metrics(m)+"_N")=n;
    S.(metrics(m)+"_CI95Low")=lo;
    S.(metrics(m)+"_CI95High")=hi;
end
S=sortrows(S,{'SensorGeometry','RequestedSNR_dB','MismatchCondition','Algorithm'});
end

function R = computeRankStability(P)
metrics=["PSNR_ROI_dB","SSIM","HD95_mm","gCNR"];
conditions=unique(P.MismatchCondition,'stable');
conditions=conditions(conditions~="BASELINE_HETERO");
geometries=unique(P.SensorGeometry,'stable');
snrs=unique(P.RequestedSNR_dB,'stable');
rows=cell(0,7);
for g=1:numel(geometries)
    for s=1:numel(snrs)
        for m=1:numel(metrics)
            metric=metrics(m);
            base=P(P.SensorGeometry==geometries(g) & P.RequestedSNR_dB==snrs(s) & ...
                P.MismatchCondition=="BASELINE_HETERO",:);
            if isempty(base), continue; end
            [algBase,valBase]=cohortMeans(base,metric);
            for c=1:numel(conditions)
                cur=P(P.SensorGeometry==geometries(g) & P.RequestedSNR_dB==snrs(s) & ...
                    P.MismatchCondition==conditions(c),:);
                [algCur,valCur]=cohortMeans(cur,metric);
                [common,ib,ic]=intersect(algBase,algCur,'stable');
                if numel(common)<3, continue; end
                descending = metric~="HD95_mm";
                rb=simpleRanks(valBase(ib),descending);
                rc=simpleRanks(valCur(ic),descending);
                rho=corr(double(rb(:)),double(rc(:)));
                reversals=countPairwiseReversals(rb,rc);
                rows(end+1,:)={geometries(g),snrs(s),metric,conditions(c),rho, ...
                    numel(common),reversals}; %#ok<AGROW>
            end
        end
    end
end
R=cell2table(rows,'VariableNames',{'SensorGeometry','RequestedSNR_dB','Metric', ...
    'MismatchCondition','SpearmanRankCorrelation','NumAlgorithms','PairwiseOrderReversals'});
if ~isempty(R)
    R.SensorGeometry=string(R.SensorGeometry); R.Metric=string(R.Metric);
    R.MismatchCondition=string(R.MismatchCondition);
end
end

function M = makePerturbationMagnitudeTable(T)
needed=["SoundSpeedMAE_mps","SoundSpeedRMSE_mps","SoundSpeedBias_mps", ...
    "DelayRMSE_vsKnownMap_ns","DelayP95Abs_vsKnownMap_ns"];
if ~all(ismember(needed,string(T.Properties.VariableNames)))
    M=table();
    return;
end
[G,condition,geometry] = findgroups(T.MismatchCondition,T.SensorGeometry);
M=table(condition,geometry,'VariableNames',{'MismatchCondition','SensorGeometry'});
for k=1:numel(needed)
    M.(needed(k))=splitapply(@meanFinite,T.(needed(k)),G);
end
M=sortrows(M,{'SensorGeometry','MismatchCondition'});
end

function files = drawDirectDeltaPSNR(S,outputDir)
direct=["DAS","CF-DAS","SCF-DAS","DMAS","DASDSF","NC-CF","NC-VDAS"];
conditions=["HOMOGENEOUS_1540","GLOBAL_MINUS_2PCT","GLOBAL_PLUS_2PCT", ...
    "MISREG_X_0P5MM","LOWRES_X4","UNDERSEGMENT_2PX"];
labels={'Homog. 1540','Global -2%','Global +2%','Shift 0.5 mm','Low-res x4','Underseg. 2 px'};
panels={"CIRCULAR_FULL",20;"CIRCULAR_FULL",0;"LIMITED_VIEW_ARC",20;"LIMITED_VIEW_ARC",0};
fig=figure('Color','w','Position',[50 50 1600 950]);
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
firstAx=[];
for p=1:4
    ax=nexttile;
    if p==1, firstAx=ax; end
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for a=1:numel(direct)
        y=nan(1,numel(conditions)); lo=y; hi=y;
        for c=1:numel(conditions)
            q=S.SensorGeometry==panels{p,1} & S.RequestedSNR_dB==panels{p,2} & ...
                S.MismatchCondition==conditions(c) & S.Algorithm==direct(a);
            if any(q)
                y(c)=S.DeltaPSNR_dB_Mean(q);
                lo(c)=S.DeltaPSNR_dB_CI95Low(q);
                hi(c)=S.DeltaPSNR_dB_CI95High(q);
            end
        end
        e1=y-lo; e2=hi-y;
        errorbar(ax,1:numel(conditions),y,e1,e2,'-o','LineWidth',1.1, ...
            'DisplayName',char(direct(a)));
    end
    yline(ax,0,':'); xlim(ax,[0.7 numel(conditions)+0.3]);
    xticks(ax,1:numel(conditions)); xticklabels(ax,labels); xtickangle(ax,25);
    ylabel(ax,'\DeltaPSNR vs known-map straight-ray baseline (dB)');
    if p==1, title(ax,'(a) Full view, 20 dB'); end
    if p==2, title(ax,'(b) Full view, 0 dB'); end
    if p==3, title(ax,'(c) Limited view, 20 dB'); end
    if p==4, title(ax,'(d) Limited view, 0 dB'); end
end
legend(firstAx,'Location','best');
title(tl,'Sensitivity of direct beamformers to sound-speed-map mismatch');
base=fullfile(outputDir,'figure_delay_mismatch_direct_delta_psnr');
exportgraphics(fig,base+".pdf",'ContentType','vector');
exportgraphics(fig,base+".png",'Resolution',240); close(fig);
files=[base+".pdf",base+".png"];
end

function files = drawAdaptiveLimitedPSNR(S,outputDir)
methods=["DAS","SM-MV","CF-SM-MV","DASDSF-SM-MV","NCVDAS-SM-MV"];
conditions=["BASELINE_HETERO","HOMOGENEOUS_1540","GLOBAL_MINUS_2PCT", ...
    "GLOBAL_PLUS_2PCT","MISREG_X_0P5MM","LOWRES_X4","UNDERSEGMENT_2PX"];
labels={'Known-map straight','Homog. 1540','Global -2%','Global +2%', ...
    'Shift 0.5 mm','Low-res x4','Underseg. 2 px'};
fig=figure('Color','w','Position',[50 50 1500 620]);
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for pp=1:2
    snr=[20 0]; ax=nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for a=1:numel(methods)
        y=nan(1,numel(conditions)); lo=y; hi=y;
        for c=1:numel(conditions)
            q=S.SensorGeometry=="LIMITED_VIEW_ARC" & S.RequestedSNR_dB==snr(pp) & ...
                S.MismatchCondition==conditions(c) & S.Algorithm==methods(a);
            if any(q)
                y(c)=S.PSNR_ROI_dB_Mean(q); lo(c)=S.PSNR_ROI_dB_CI95Low(q);
                hi(c)=S.PSNR_ROI_dB_CI95High(q);
            end
        end
        errorbar(ax,1:numel(conditions),y,y-lo,hi-y,'-o','LineWidth',1.1, ...
            'DisplayName',char(methods(a)));
    end
    xlim(ax,[0.7 numel(conditions)+0.3]); xticks(ax,1:numel(conditions));
    xticklabels(ax,labels); xtickangle(ax,28); ylabel(ax,'PSNR (dB)');
    title(ax,sprintf('Limited view, %d dB',snr(pp)));
    if pp==1, legend(ax,'Location','best'); end
end
title(tl,'Adaptive-beamforming comparison under sound-speed-map mismatch');
base=fullfile(outputDir,'figure_delay_mismatch_adaptive_limited_psnr');
exportgraphics(fig,base+".pdf",'ContentType','vector');
exportgraphics(fig,base+".png",'Resolution',240); close(fig);
files=[base+".pdf",base+".png"];
end

function [algs,vals]=cohortMeans(T,metric)
algs=unique(T.Algorithm,'stable'); vals=nan(numel(algs),1);
for k=1:numel(algs), vals(k)=mean(T.(metric)(T.Algorithm==algs(k)),'omitnan'); end
end
function r=simpleRanks(v,descending)
v=double(v(:)); if descending, [~,o]=sort(v,'descend'); else, [~,o]=sort(v,'ascend'); end
r=zeros(size(v)); r(o)=1:numel(v);
end
function n=countPairwiseReversals(a,b)
n=0; for i=1:numel(a)-1, for j=i+1:numel(a), if sign(a(i)-a(j))*sign(b(i)-b(j))<0, n=n+1; end, end, end
end
function ci=bootstrapPair(v,B)
v=double(v(:)); v=v(isfinite(v));
if isempty(v), ci=[NaN NaN]; return; end
if numel(v)==1, ci=[v v]; return; end
n=numel(v); idx=randi(n,n,B); bm=mean(v(idx),1); ci=[percentileFinite(bm,2.5),percentileFinite(bm,97.5)];
end
function m=meanFinite(v), v=double(v(:)); v=v(isfinite(v)); if isempty(v), m=NaN; else, m=mean(v); end, end
function s=stdFinite(v), v=double(v(:)); v=v(isfinite(v)); if numel(v)<2, s=NaN; else, s=std(v,0); end, end
function n=countFinite(v), n=nnz(isfinite(v)); end
function q=percentileFinite(x,pct)
x=sort(double(x(isfinite(x)))); if isempty(x), q=NaN; return; end
pos=1+(numel(x)-1)*pct/100; lo=floor(pos); hi=ceil(pos);
if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end
