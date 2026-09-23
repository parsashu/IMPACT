function outputs = analyze_adaptive_cohort( ...
    newAdaptiveCsv, mainResultsCsv, outputDir)
% analyze_adaptive_cohort
% -------------------------------------------------------------------------
% Combine the original five-phantom adaptive cohort with the six new v9.6
% adaptive-only rows, then compare each adaptive method with DAS across all
% 11 phantoms in LIMITED_VIEW_ARC.
%
% Statistical unit:
%   PHANTOM. Repeated noise realizations are averaged within each phantom
%   before effect-size or hypothesis calculations.
%
% Reported quantities:
%   - paired raw mean difference (candidate - DAS)
%   - percentile 95% phantom bootstrap CI of the paired mean difference
%   - paired standardized effect size dz
%   - median and IQR of paired differences
%   - number of phantoms with positive / negative difference
%   - exact sign-flip permutation p-value for n <= 16
%   - Benjamini-Hochberg FDR q-value
%
% No equivalence/non-inferiority claim is made. A non-significant test is
% explicitly NOT interpreted as equivalence.
% -------------------------------------------------------------------------

if nargin<1 || isempty(newAdaptiveCsv)
    newAdaptiveCsv = fullfile(pwd,'ADAPTIVE_11PHANTOM_COMPLETION', ...
        'adaptive_11phantom_results.csv');
end
if nargin<2 || isempty(mainResultsCsv)
    mainResultsCsv = fullfile(pwd,'benchmark_v93_reviewer_revision', ...
        'all_reconstruction_results_PARTIAL.csv');
end
if nargin<3 || isempty(outputDir)
    outputDir = fullfile(pwd,'ADAPTIVE_11PHANTOM_ANALYSIS');
end
if ~isfolder(outputDir), mkdir(outputDir); end

if ~isfile(newAdaptiveCsv)
    error('PACT:MissingAdaptiveCompletion','Missing: %s',newAdaptiveCsv);
end
if ~isfile(mainResultsCsv)
    error('PACT:MissingMainResults','Missing: %s',mainResultsCsv);
end

N = readtable(newAdaptiveCsv,'TextType','string');
M = readtable(mainResultsCsv,'TextType','string');

geometry = "LIMITED_VIEW_ARC";
adaptiveAlgorithms = ["SM-MV","CF-SM-MV","DASDSF-SM-MV","NCVDAS-SM-MV"];
snrValues = [Inf,20,10,0,-5];
allPhantoms = ["Circle","Ring","Two_Circles","Ellipse","Rectangle","Triangle", ...
    "Shepp_Logan","Complex","Medical","Letter_A","Vascular_Tree"];
bootstrapCount = 4000;
statisticsSeed = 661177;

% Original adaptive rows: retain only the same limited-view schedule.
O = M(M.SensorGeometry==geometry & ...
    ismember(M.Algorithm,adaptiveAlgorithms),:);

% Convert the new table to the common metric schema used below.
commonNames = ["Phantom","SensorGeometry","Algorithm","RequestedSNR_dB", ...
    "NoiseRealization","PSNR_ROI_dB","SSIM","Dice","HD95_mm","gCNR", ...
    "ArtifactToTarget_dB"];
O = O(:,commonNames);
N = N(:,commonNames);

A = [O;N];
A = uniqueAdaptiveRowsForAnalysis(A);

% Keep exactly the schedule used for the adaptive cohort.
keep = ismember(A.Phantom,allPhantoms) & A.SensorGeometry==geometry & ...
    ismember(A.Algorithm,adaptiveAlgorithms);
A = A(keep,:);

% DAS references from the original 11-phantom main benchmark.
D = M(M.SensorGeometry==geometry & M.Algorithm=="DAS" & ...
    ismember(M.Phantom,allPhantoms),commonNames);
D = uniqueAdaptiveRowsForAnalysis(D);

% Strict completeness checks before statistics.
for phantom = allPhantoms
    for snrDb = snrValues
        reps = 5; if isinf(snrDb), reps=1; end
        for realization = 1:reps
            idxD = D.Phantom==phantom & snrVectorEqualLocal(D.RequestedSNR_dB,snrDb) & ...
                D.NoiseRealization==realization;
            if nnz(idxD)~=1
                error('PACT:IncompleteDAS', ...
                    'DAS completeness failure: %s SNR=%s rep=%d.', ...
                    phantom,formatSNRLocal(snrDb),realization);
            end
            for algorithm = adaptiveAlgorithms
                idxA = A.Phantom==phantom & A.Algorithm==algorithm & ...
                    snrVectorEqualLocal(A.RequestedSNR_dB,snrDb) & ...
                    A.NoiseRealization==realization;
                if nnz(idxA)~=1
                    error('PACT:IncompleteAdaptive', ...
                        'Adaptive completeness failure: %s | %s | SNR=%s | rep=%d.', ...
                        phantom,algorithm,formatSNRLocal(snrDb),realization);
                end
            end
        end
    end
end

expectedAdaptiveRows = numel(allPhantoms)*(1+4*5)*numel(adaptiveAlgorithms);
if height(A)~=expectedAdaptiveRows
    error('PACT:AdaptiveRowCount', ...
        'Expected %d combined adaptive rows; found %d.', ...
        expectedAdaptiveRows,height(A));
end
fprintf('Combined adaptive cohort completeness: PASS (%d rows).\n',height(A));

writetable(A,fullfile(outputDir,'adaptive_11phantom_combined_rows.csv'));

%% Phantom-level means
metrics = ["PSNR_ROI_dB","SSIM","Dice","HD95_mm","gCNR","ArtifactToTarget_dB"];
phantomRows = cell(0,11);
for phantom = allPhantoms
    for snrDb = snrValues
        for algorithm = ["DAS",adaptiveAlgorithms]
            if algorithm=="DAS"
                T = D;
            else
                T = A;
            end
            idx = T.Phantom==phantom & T.Algorithm==algorithm & ...
                snrVectorEqualLocal(T.RequestedSNR_dB,snrDb);
            if ~any(idx), continue; end
            values = cell(1,numel(metrics));
            for m = 1:numel(metrics)
                values{m} = mean(T.(metrics(m))(idx),'omitnan');
            end
            phantomRows(end+1,:) = [{phantom,geometry,algorithm,snrDb},values]; %#ok<AGROW>
        end
    end
end

P = cell2table(phantomRows,'VariableNames', ...
    [{'Phantom','SensorGeometry','Algorithm','RequestedSNR_dB'},cellstr(metrics)]);
P.Phantom = string(P.Phantom);
P.SensorGeometry = string(P.SensorGeometry);
P.Algorithm = string(P.Algorithm);
writetable(P,fullfile(outputDir,'adaptive_11phantom_phantom_means.csv'));

%% Paired effect sizes and exact tests
higherIsBetter = [true,true,true,false,true,false];
rng(statisticsSeed,'twister');
records = cell(0,1);
row = 1;

for snrDb = snrValues
    baseline = P(P.Algorithm=="DAS" & ...
        snrVectorEqualLocal(P.RequestedSNR_dB,snrDb),:);

    for algorithm = adaptiveAlgorithms
        candidate = P(P.Algorithm==algorithm & ...
            snrVectorEqualLocal(P.RequestedSNR_dB,snrDb),:);

        [commonPhantoms,ib,ic] = intersect(baseline.Phantom,candidate.Phantom,'stable');
        if numel(commonPhantoms)~=11
            error('PACT:PairingFailure', ...
                'Expected 11 paired phantoms for %s at SNR=%s; found %d.', ...
                algorithm,formatSNRLocal(snrDb),numel(commonPhantoms));
        end

        for m = 1:numel(metrics)
            b = baseline.(metrics(m))(ib);
            c = candidate.(metrics(m))(ic);
            rawDiff = c-b;
            if higherIsBetter(m)
                orientedDiff = rawDiff;
            else
                orientedDiff = -rawDiff;
            end

            [ciLow,ciHigh] = bootstrapMeanCILocal(rawDiff,bootstrapCount);
            p = exactSignFlipPLocal(orientedDiff);
            sd = std(orientedDiff,0,'omitnan');
            if sd<=eps
                dz = NaN;
            else
                dz = mean(orientedDiff,'omitnan')/sd;
            end

            rec = struct();
            rec.EvaluationUnit = "PHANTOM";
            rec.SensorGeometry = geometry;
            rec.RequestedSNR_dB = snrDb;
            rec.Algorithm = algorithm;
            rec.Baseline = "DAS";
            rec.Metric = metrics(m);
            rec.NumPhantoms = numel(commonPhantoms);
            rec.MeanCandidate = mean(c,'omitnan');
            rec.MeanBaseline = mean(b,'omitnan');
            rec.MeanDifference_CandidateMinusDAS = mean(rawDiff,'omitnan');
            rec.CI95Low = ciLow;
            rec.CI95High = ciHigh;
            rec.MedianDifference = median(rawDiff,'omitnan');
            rec.IQRDifference = percentileLinearLocal(rawDiff,75)- ...
                percentileLinearLocal(rawDiff,25);
            rec.NumPositiveRawDifference = nnz(rawDiff>0);
            rec.NumNegativeRawDifference = nnz(rawDiff<0);
            rec.OrientedEffectSizeDz = dz;
            rec.ExactSignFlipP = p;
            rec.HigherIsBetter = higherIsBetter(m);
            records{row,1} = rec; %#ok<AGROW>
            row = row+1;
        end
    end
end

E = struct2table(vertcat(records{:}));
E.FDR_Q = benjaminiHochbergLocal(E.ExactSignFlipP);
writetable(E,fullfile(outputDir,'adaptive_11phantom_effect_sizes_vs_DAS.csv'));

% Compact manuscript table: structural metrics most relevant to the reviewer.
keyMetrics = ["PSNR_ROI_dB","SSIM","Dice","HD95_mm"];
K = E(ismember(E.Metric,keyMetrics),:);
writetable(K,fullfile(outputDir,'Table_adaptive_effects_vs_DAS.csv'));

%% Absolute 11-phantom summary
summaryRecords = cell(0,1);
row = 1;
for snrDb = snrValues
    for algorithm = ["DAS",adaptiveAlgorithms]
        B = P(P.Algorithm==algorithm & ...
            snrVectorEqualLocal(P.RequestedSNR_dB,snrDb),:);
        for m = 1:numel(metrics)
            v = B.(metrics(m));
            rec = struct();
            rec.SensorGeometry = geometry;
            rec.RequestedSNR_dB = snrDb;
            rec.Algorithm = algorithm;
            rec.Metric = metrics(m);
            rec.NumPhantoms = nnz(isfinite(v));
            rec.Mean = mean(v,'omitnan');
            rec.Std = std(v,0,'omitnan');
            rec.Median = median(v,'omitnan');
            rec.IQR = percentileLinearLocal(v,75)-percentileLinearLocal(v,25);
            [rec.CI95Low,rec.CI95High] = bootstrapMeanCILocal(v,bootstrapCount);
            summaryRecords{row,1} = rec; %#ok<AGROW>
            row = row+1;
        end
    end
end
S = struct2table(vertcat(summaryRecords{:}));
writetable(S,fullfile(outputDir,'adaptive_11phantom_absolute_summary.csv'));

%% Figure: raw paired effect with 95% phantom bootstrap CI
fig = figure('Color','w','Position',[50 40 1550 1050]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

plotEffectPanel(nexttile(tl,1),E,snrValues,adaptiveAlgorithms, ...
    "PSNR_ROI_dB","\DeltaPSNR [dB]","(a) PSNR: adaptive minus DAS");
plotEffectPanel(nexttile(tl,2),E,snrValues,adaptiveAlgorithms, ...
    "SSIM","\DeltaSSIM","(b) SSIM: adaptive minus DAS");
plotEffectPanel(nexttile(tl,3),E,snrValues,adaptiveAlgorithms, ...
    "Dice","\DeltaDice","(c) Dice: adaptive minus DAS");
plotEffectPanel(nexttile(tl,4),E,snrValues,adaptiveAlgorithms, ...
    "gCNR","\DeltagCNR","(d) gCNR: adaptive minus DAS");

title(tl, ...
    "Limited-view adaptive beamforming: paired effects across all 11 phantoms", ...
    'FontWeight','bold','FontSize',18);

pngFile = fullfile(outputDir,'Figure_adaptive_11phantom_effects.png');
pdfFile = fullfile(outputDir,'Figure_adaptive_11phantom_effects.pdf');
exportgraphics(fig,pngFile,'Resolution',300);
exportgraphics(fig,pdfFile,'ContentType','vector');
close(fig);

%% A second diagnostic figure with all individual phantom differences for SSIM
fig2 = figure('Color','w','Position',[70 70 1300 650]);
ax = axes(fig2); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plotPhantomSSIMDifferences(ax,P,snrValues,adaptiveAlgorithms);
title(ax,'Phantom-level SSIM differences relative to DAS (limited view)');
ylabel(ax,'\DeltaSSIM');
xlabel(ax,'Input SNR [dB]');
pngFile2 = fullfile(outputDir,'Figure_adaptive_11phantom_phantom_SSIM.png');
pdfFile2 = fullfile(outputDir,'Figure_adaptive_11phantom_phantom_SSIM.pdf');
exportgraphics(fig2,pngFile2,'Resolution',300);
exportgraphics(fig2,pdfFile2,'ContentType','vector');
close(fig2);

%% Console report
fprintf('\n===============================================================\n');
fprintf(' 11-PHANTOM ADAPTIVE EFFECT-SIZE ANALYSIS COMPLETE\n');
fprintf('===============================================================\n');
for snrDb = snrValues
    fprintf('\nSNR %s dB\n',formatSNRLocal(snrDb));
    for algorithm = adaptiveAlgorithms
        rS = E(E.Algorithm==algorithm & E.Metric=="SSIM" & ...
            snrVectorEqualLocal(E.RequestedSNR_dB,snrDb),:);
        rP = E(E.Algorithm==algorithm & E.Metric=="PSNR_ROI_dB" & ...
            snrVectorEqualLocal(E.RequestedSNR_dB,snrDb),:);
        fprintf(['  %-16s | dSSIM=%+.4f [%+.4f,%+.4f], q=%.4g | ' ...
            'dPSNR=%+.3f [%+.3f,%+.3f] dB, q=%.4g\n'], ...
            algorithm,rS.MeanDifference_CandidateMinusDAS, ...
            rS.CI95Low,rS.CI95High,rS.FDR_Q, ...
            rP.MeanDifference_CandidateMinusDAS, ...
            rP.CI95Low,rP.CI95High,rP.FDR_Q);
    end
end

outputs = struct();
outputs.combinedAdaptiveRows = A;
outputs.phantomMeans = P;
outputs.effectSizes = E;
outputs.absoluteSummary = S;
outputs.figurePNG = string(pngFile);
outputs.figurePDF = string(pdfFile);
outputs.phantomSSIMFigurePNG = string(pngFile2);
outputs.phantomSSIMFigurePDF = string(pdfFile2);
outputs.outputDir = string(outputDir);
end

%% ========================================================================
function T = uniqueAdaptiveRowsForAnalysis(T)
if isempty(T), return; end
keys = T.Phantom+"|"+T.SensorGeometry+"|"+T.Algorithm+"|"+ ...
    string(T.RequestedSNR_dB)+"|"+string(T.NoiseRealization);
[~,idx] = unique(keys,'stable');
T = T(idx,:);
end

%% ========================================================================
function mask = snrVectorEqualLocal(values,target)
if isinf(target)
    mask = isinf(values) & sign(values)==sign(target);
else
    mask = abs(values-target)<1e-9;
end
end

%% ========================================================================
function text = formatSNRLocal(value)
if isinf(value)
    if value>0, text='Inf'; else, text='-Inf'; end
else
    text=sprintf('%.0f',value);
end
end

%% ========================================================================
function [low,high] = bootstrapMeanCILocal(values,count)
v = values(:);
v = v(isfinite(v));
n = numel(v);
if n<2
    low=NaN; high=NaN; return;
end
means = zeros(count,1);
chunk = 500;
filled = 0;
while filled<count
    current = min(chunk,count-filled);
    idx = randi(n,n,current);
    means(filled+(1:current)) = mean(v(idx),1).';
    filled = filled+current;
end
low = percentileLinearLocal(means,2.5);
high = percentileLinearLocal(means,97.5);
end

%% ========================================================================
function p = exactSignFlipPLocal(orientedDifferences)
d = orientedDifferences(:);
d = d(isfinite(d));
n = numel(d);
if n==0, p=NaN; return; end
if n>16
    error('Exact sign-flip helper is configured for n<=16; got n=%d.',n);
end
observed = abs(mean(d));
numberPatterns = 2^n;
integers = uint32((0:numberPatterns-1)');
signs = ones(n,numberPatterns);
for bitIndex = 1:n
    bits = bitget(integers,bitIndex).';
    signs(bitIndex,:) = 2*double(bits)-1;
end
permutedMeans = abs(mean(d.*signs,1));
p = nnz(permutedMeans>=observed-10*eps)/numberPatterns;
end

%% ========================================================================
function q = benjaminiHochbergLocal(p)
p = p(:);
q = NaN(size(p));
valid = find(isfinite(p));
if isempty(valid), return; end
[pSorted,order] = sort(p(valid));
m = numel(pSorted);
qSorted = pSorted.*m./(1:m)';
for k=m-1:-1:1
    qSorted(k)=min(qSorted(k),qSorted(k+1));
end
qSorted=min(qSorted,1);
temp=zeros(m,1);
temp(order)=qSorted;
q(valid)=temp;
end

%% ========================================================================
function value = percentileLinearLocal(x,percent)
x = sort(x(isfinite(x)));
if isempty(x), value=NaN; return; end
if numel(x)==1, value=x; return; end
position = 1 + (numel(x)-1)*percent/100;
lo=floor(position); hi=ceil(position);
if lo==hi
    value=x(lo);
else
    value=x(lo)+(position-lo)*(x(hi)-x(lo));
end
end

%% ========================================================================
function plotEffectPanel(ax,E,snrValues,algorithms,metric,ylab,ttl)
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
yline(ax,0,'--','DAS reference','HandleVisibility','off');

x = 1:numel(snrValues);
offsets = linspace(-0.24,0.24,numel(algorithms));
markers = {'o','s','^','d'};

for a=1:numel(algorithms)
    means=nan(size(x)); lo=nan(size(x)); hi=nan(size(x));
    for s=1:numel(snrValues)
        idx = E.Algorithm==algorithms(a) & E.Metric==metric & ...
            snrVectorEqualLocal(E.RequestedSNR_dB,snrValues(s));
        if ~any(idx), continue; end
        means(s)=E.MeanDifference_CandidateMinusDAS(idx);
        lo(s)=means(s)-E.CI95Low(idx);
        hi(s)=E.CI95High(idx)-means(s);
    end
    errorbar(ax,x+offsets(a),means,lo,hi,'-','LineWidth',1.25, ...
        'Marker',markers{a},'MarkerSize',6,'DisplayName',algorithms(a));
end
xticks(ax,x);
labels=cell(size(snrValues));
for k=1:numel(snrValues), labels{k}=formatSNRLocal(snrValues(k)); end
xticklabels(ax,labels);
xlabel(ax,'Input SNR [dB]');
ylabel(ax,ylab);
title(ax,ttl,'FontWeight','bold');
legend(ax,'Location','best','Box','off');
set(ax,'FontSize',11,'LineWidth',0.9);
end

%% ========================================================================
function plotPhantomSSIMDifferences(ax,P,snrValues,algorithms)
colors = lines(numel(algorithms));
x = 1:numel(snrValues);
for a=1:numel(algorithms)
    for s=1:numel(snrValues)
        b=P(P.Algorithm=="DAS" & ...
            snrVectorEqualLocal(P.RequestedSNR_dB,snrValues(s)),:);
        c=P(P.Algorithm==algorithms(a) & ...
            snrVectorEqualLocal(P.RequestedSNR_dB,snrValues(s)),:);
        [~,ib,ic]=intersect(b.Phantom,c.Phantom,'stable');
        diff=c.SSIM(ic)-b.SSIM(ib);
        xx=x(s)+linspace(-0.08,0.08,numel(diff));
        scatter(ax,xx,diff,20,colors(a,:),'filled', ...
            'MarkerFaceAlpha',0.45,'HandleVisibility','off');
        plot(ax,x(s),mean(diff),'Marker','o','MarkerSize',8, ...
            'Color',colors(a,:),'MarkerFaceColor',colors(a,:), ...
            'HandleVisibility','off');
    end
    plot(ax,nan,nan,'-o','Color',colors(a,:),'DisplayName',algorithms(a));
end
yline(ax,0,'--','DAS reference','HandleVisibility','off');
xticks(ax,x);
labels=cell(size(snrValues));
for k=1:numel(snrValues), labels{k}=formatSNRLocal(snrValues(k)); end
xticklabels(ax,labels);
legend(ax,'Location','best','Box','off');
end
