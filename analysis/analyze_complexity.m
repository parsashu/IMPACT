
function outputs = analyze_complexity( ...
    resolutionFile, complexityFile, slopeFile, outputDir)
% PACT_PLOT_RESOLUTION_COMPLEXITY
%
% Reads the final PACT CSV files and generates two manuscript-ready figures:
%
%   Figure 1: Point-target FWHM and peak-sidelobe level
%             (full-view and limited-view; 20 and 0 dB)
%
%   Figure 2: Direct SS-MV versus sequential SM-MV computational scaling
%             (runtime, fitted slopes, and SM/direct runtime ratio)
%
% No simulation or reconstruction is repeated.
%
% Required files:
%   resolution_results.csv
%   mv_complexity_results.csv
%   mv_complexity_slopes.csv
%
% Example:
%   outputs = analyze_complexity( ...
%       "resolution_results.csv", ...
%       "mv_complexity_results.csv", ...
%       "mv_complexity_slopes.csv", ...
%       "ARTICLE_FIGURES");
%
% Outputs:
%   ARTICLE_FIGURES/resolution_summary.csv
%   ARTICLE_FIGURES/mv_runtime_ratio.csv
%   ARTICLE_FIGURES/figure_resolution.pdf/png/fig
%   ARTICLE_FIGURES/figure_mv_complexity.pdf/png/fig
%   ARTICLE_FIGURES/latex_include_snippet.tex
%
% MATLAB R2025a compatible.

    if nargin < 1 || strlength(string(resolutionFile)) == 0
        resolutionFile = "resolution_results.csv";
    end
    if nargin < 2 || strlength(string(complexityFile)) == 0
        complexityFile = "mv_complexity_results.csv";
    end
    if nargin < 3 || strlength(string(slopeFile)) == 0
        slopeFile = "mv_complexity_slopes.csv";
    end
    if nargin < 4 || strlength(string(outputDir)) == 0
        outputDir = "PACT_ARTICLE_FIGURES";
    end

    resolutionFile = string(resolutionFile);
    complexityFile = string(complexityFile);
    slopeFile = string(slopeFile);
    outputDir = string(outputDir);

    assertFile(resolutionFile);
    assertFile(complexityFile);
    assertFile(slopeFile);
    ensureFolder(outputDir);

    fprintf('\n============================================================\n');
    fprintf(' PACT RESOLUTION AND MV COMPLEXITY FIGURES\n');
    fprintf(' Post-processing only; no reconstruction is repeated\n');
    fprintf('============================================================\n');

    resolutionRaw = readtable(resolutionFile, ...
        'VariableNamingRule','preserve');
    complexityRaw = readtable(complexityFile, ...
        'VariableNamingRule','preserve');
    slopesRaw = readtable(slopeFile, ...
        'VariableNamingRule','preserve');

    validateResolutionColumns(resolutionRaw);
    validateComplexityColumns(complexityRaw);
    validateSlopeColumns(slopesRaw);

    resolutionRaw.SensorGeometry = string(resolutionRaw.SensorGeometry);
    resolutionRaw.Algorithm = string(resolutionRaw.Algorithm);
    complexityRaw.Method = string(complexityRaw.Method);
    slopesRaw.Method = string(slopesRaw.Method);

    resolutionSummary = summarizeResolution(resolutionRaw);
    writetable(resolutionSummary, ...
        fullfile(outputDir,'resolution_summary.csv'));

    runtimeRatio = makeRuntimeRatioTable(complexityRaw);
    writetable(runtimeRatio, ...
        fullfile(outputDir,'mv_runtime_ratio.csv'));

    resolutionFiles = drawResolutionFigure( ...
        resolutionSummary, outputDir);

    complexityFiles = drawComplexityFigure( ...
        complexityRaw, slopesRaw, runtimeRatio, outputDir);

    latexFile = writeLatexSnippet(outputDir);

    outputs = struct;
    outputs.resolutionFile = resolutionFile;
    outputs.complexityFile = complexityFile;
    outputs.slopeFile = slopeFile;
    outputs.outputDir = outputDir;
    outputs.resolutionSummary = resolutionSummary;
    outputs.runtimeRatio = runtimeRatio;
    outputs.resolutionFigureFiles = resolutionFiles;
    outputs.complexityFigureFiles = complexityFiles;
    outputs.latexSnippet = string(latexFile);

    save(fullfile(outputDir,'plot_outputs.mat'), ...
        'outputs','-v7.3');

    assignin('base','PACT_resolutionSummary',resolutionSummary);
    assignin('base','PACT_mvRuntimeRatio',runtimeRatio);
    assignin('base','PACT_plotOutputs',outputs);

    fprintf('\nCompleted.\nOutput folder:\n  %s\n',outputDir);
end

% ========================================================================
function S = summarizeResolution(T)
% Mean and SD over point locations and noise realizations.

    [G,geometry,snr,algorithm] = findgroups( ...
        string(T.SensorGeometry), ...
        T.RequestedSNR_dB, ...
        string(T.Algorithm));

    S = table(geometry,snr,algorithm, ...
        'VariableNames',{'SensorGeometry','RequestedSNR_dB','Algorithm'});

    S.N = splitapply(@numel,T.FWHM_Radial_mm,G);

    S.FWHM_Radial_Mean_mm = ...
        splitapply(@meanFinite,T.FWHM_Radial_mm,G);
    S.FWHM_Radial_SD_mm = ...
        splitapply(@stdFinite,T.FWHM_Radial_mm,G);

    S.FWHM_Tangential_Mean_mm = ...
        splitapply(@meanFinite,T.FWHM_Tangential_mm,G);
    S.FWHM_Tangential_SD_mm = ...
        splitapply(@stdFinite,T.FWHM_Tangential_mm,G);

    S.PSL_Radial_Mean_dB = ...
        splitapply(@meanFinite,T.PSL_Radial_dB,G);
    S.PSL_Radial_SD_dB = ...
        splitapply(@stdFinite,T.PSL_Radial_dB,G);

    S.PSL_Tangential_Mean_dB = ...
        splitapply(@meanFinite,T.PSL_Tangential_dB,G);
    S.PSL_Tangential_SD_dB = ...
        splitapply(@stdFinite,T.PSL_Tangential_dB,G);

    S.Localization_Mean_mm = ...
        splitapply(@meanFinite,T.PeakLocalizationError_mm,G);
    S.Localization_SD_mm = ...
        splitapply(@stdFinite,T.PeakLocalizationError_mm,G);

    S = sortrows(S, ...
        {'SensorGeometry','RequestedSNR_dB','Algorithm'});
end

% ========================================================================
function R = makeRuntimeRatioTable(T)

    direct = T(T.Method=="SS-MV Direct",:);
    sm = T(T.Method=="SM-MV",:);

    direct = direct(:,{ ...
        'NumSensors_M','SubarrayLength_L', ...
        'MedianRuntime_Sec','RuntimeIQR_Sec'});

    sm = sm(:,{ ...
        'NumSensors_M','SubarrayLength_L', ...
        'MedianRuntime_Sec','RuntimeIQR_Sec'});

    direct.Properties.VariableNames = { ...
        'NumSensors_M','SubarrayLength_L', ...
        'DirectMedianRuntime_Sec','DirectRuntimeIQR_Sec'};

    sm.Properties.VariableNames = { ...
        'NumSensors_M','SubarrayLength_L', ...
        'SMMedianRuntime_Sec','SMRuntimeIQR_Sec'};

    R = innerjoin(direct,sm, ...
        'Keys',{'NumSensors_M','SubarrayLength_L'});

    R.DirectMedianRuntime_ms = 1000*R.DirectMedianRuntime_Sec;
    R.SMMedianRuntime_ms = 1000*R.SMMedianRuntime_Sec;

    R.DirectRuntimeIQR_ms = 1000*R.DirectRuntimeIQR_Sec;
    R.SMRuntimeIQR_ms = 1000*R.SMRuntimeIQR_Sec;

    R.SM_to_Direct_Ratio = ...
        R.SMMedianRuntime_Sec ./ R.DirectMedianRuntime_Sec;

    R = sortrows(R,{'SubarrayLength_L','NumSensors_M'});
end

% ========================================================================
function files = drawResolutionFigure(S,outputDir)

    fullMethods = ["DAS","CF-DAS","SCF-DAS","DMAS","DASDSF"];

    limitedMethods = [ ...
        "DAS","DMAS","DASDSF","SM-MV", ...
        "CF-SM-MV","DASDSF-SM-MV","NCVDAS-SM-MV"];

    snrValues = [20 0];

    fig = figure( ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[60 40 1900 1120]);

    tl = tiledlayout(fig,2,2, ...
        'TileSpacing','loose', ...
        'Padding','loose');

    ax1 = nexttile(tl,1);
    legendHandles = plotResolutionPanel( ...
        ax1,S,"CIRCULAR_FULL",fullMethods,snrValues,"FWHM");
    title(ax1,'(a) Full-view FWHM', ...
        'FontWeight','bold');

    ax2 = nexttile(tl,2);
    plotResolutionPanel( ...
        ax2,S,"LIMITED_VIEW_ARC",limitedMethods,snrValues,"FWHM");
    title(ax2,'(b) Limited-view FWHM', ...
        'FontWeight','bold');

    ax3 = nexttile(tl,3);
    plotResolutionPanel( ...
        ax3,S,"CIRCULAR_FULL",fullMethods,snrValues,"PSL");
    title(ax3,'(c) Full-view peak-sidelobe level', ...
        'FontWeight','bold');

    ax4 = nexttile(tl,4);
    plotResolutionPanel( ...
        ax4,S,"LIMITED_VIEW_ARC",limitedMethods,snrValues,"PSL");
    title(ax4,'(d) Limited-view peak-sidelobe level', ...
        'FontWeight','bold');

    % Place the legend on the first axes only. This is compatible with
    % MATLAB releases that do not support legends attached to tiledlayout.
    lgd = legend(ax1,legendHandles, ...
        {'Radial, 20 dB','Tangential, 20 dB', ...
         'Radial, 0 dB','Tangential, 0 dB'}, ...
        'Location','northoutside', ...
        'Orientation','horizontal', ...
        'NumColumns',4, ...
        'Box','off', ...
        'FontName','Times New Roman', ...
        'FontSize',10);


    title(tl, ...
        'Point-target resolution and sidelobe behavior', ...
        'FontName','Times New Roman', ...
        'FontSize',19, ...
        'FontWeight','bold');

    files = exportFigureSet( ...
        fig,outputDir,'figure_resolution');

    close(fig);
end

% ========================================================================
function handles = plotResolutionPanel( ...
        ax,S,geometry,methods,snrValues,metricType)

    hold(ax,'on');
    box(ax,'on');
    grid(ax,'on');

    x = 1:numel(methods);
    offsets = [-0.21 -0.07 0.07 0.21];
    colors = lines(4);
    markers = {'o','s','o','s'};
    filled = [true true false false];

    if metricType=="FWHM"
        meanColumns = [ ...
            "FWHM_Radial_Mean_mm", ...
            "FWHM_Tangential_Mean_mm", ...
            "FWHM_Radial_Mean_mm", ...
            "FWHM_Tangential_Mean_mm"];

        sdColumns = [ ...
            "FWHM_Radial_SD_mm", ...
            "FWHM_Tangential_SD_mm", ...
            "FWHM_Radial_SD_mm", ...
            "FWHM_Tangential_SD_mm"];

        ylabel(ax,'FWHM (mm)', ...
            'FontWeight','bold');
    else
        meanColumns = [ ...
            "PSL_Radial_Mean_dB", ...
            "PSL_Tangential_Mean_dB", ...
            "PSL_Radial_Mean_dB", ...
            "PSL_Tangential_Mean_dB"];

        sdColumns = [ ...
            "PSL_Radial_SD_dB", ...
            "PSL_Tangential_SD_dB", ...
            "PSL_Radial_SD_dB", ...
            "PSL_Tangential_SD_dB"];

        ylabel(ax,'Peak-sidelobe level (dB)', ...
            'FontWeight','bold');
        yline(ax,0,':','Color',[0.35 0.35 0.35]);
    end

    seriesSNR = [snrValues(1),snrValues(1), ...
                 snrValues(2),snrValues(2)];

    handles = gobjects(4,1);

    for k = 1:4
        [y,e] = extractResolutionSeries( ...
            S,geometry,seriesSNR(k),methods, ...
            meanColumns(k),sdColumns(k));

        if filled(k)
            markerFace = colors(k,:);
        else
            markerFace = 'w';
        end

        handles(k) = errorbar( ...
            ax,x+offsets(k),y,e, ...
            'LineStyle','none', ...
            'Marker',markers{k}, ...
            'MarkerSize',7, ...
            'MarkerFaceColor',markerFace, ...
            'MarkerEdgeColor',colors(k,:), ...
            'Color',colors(k,:), ...
            'LineWidth',1.15, ...
            'CapSize',7);
    end

    ax.XLim = [0.45 numel(methods)+0.55];
    ax.XTick = x;
    ax.XTickLabel = cellstr(methods);
    ax.XTickLabelRotation = 31;

    ax.FontName = 'Times New Roman';
    ax.FontSize = 11;
    ax.LineWidth = 0.9;
    ax.GridAlpha = 0.20;

    if metricType=="FWHM"
        allValues = [];
        for k = 1:4
            [y,e] = extractResolutionSeries( ...
                S,geometry,seriesSNR(k),methods, ...
                meanColumns(k),sdColumns(k));
            allValues = [allValues; y(:)+e(:)]; %#ok<AGROW>
        end
        upper = max(allValues(isfinite(allValues)));
        if isempty(upper)
            upper = 1;
        end
        ylim(ax,[0 1.12*upper]);
    else
        allLow = [];
        allHigh = [];
        for k = 1:4
            [y,e] = extractResolutionSeries( ...
                S,geometry,seriesSNR(k),methods, ...
                meanColumns(k),sdColumns(k));
            allLow = [allLow; y(:)-e(:)]; %#ok<AGROW>
            allHigh = [allHigh; y(:)+e(:)]; %#ok<AGROW>
        end

        low = min(allLow(isfinite(allLow)));
        high = max(allHigh(isfinite(allHigh)));

        if isempty(low), low = -15; end
        if isempty(high), high = 0; end

        span = max(high-low,1);
        ylim(ax,[low-0.08*span, min(1,high+0.12*span)]);
    end
end

% ========================================================================
function [y,e] = extractResolutionSeries( ...
        S,geometry,snrValue,methods,meanColumn,sdColumn)

    y = NaN(1,numel(methods));
    e = NaN(1,numel(methods));

    for m = 1:numel(methods)
        mask = ...
            S.SensorGeometry==geometry & ...
            abs(S.RequestedSNR_dB-snrValue)<1e-12 & ...
            S.Algorithm==methods(m);

        if any(mask)
            row = find(mask,1,'first');
            y(m) = S.(meanColumn)(row);
            e(m) = S.(sdColumn)(row);
        end
    end
end

% ========================================================================
function files = drawComplexityFigure( ...
        T,slopes,R,outputDir)

    subarrayLengths = unique(T.SubarrayLength_L,'sorted')';

    fig = figure( ...
        'Color','w', ...
        'Units','pixels', ...
        'Position',[70 70 1850 660]);

    tl = tiledlayout(fig,1,3, ...
        'TileSpacing','loose', ...
        'Padding','loose');

    % ------------------------------------------------------------
    % Panel (a): runtime scaling
    % ------------------------------------------------------------
    ax1 = nexttile(tl,1);
    hold(ax1,'on');
    box(ax1,'on');
    grid(ax1,'on');

    colors = lines(numel(subarrayLengths));
    runtimeHandles = gobjects(2*numel(subarrayLengths),1);
    runtimeLabels = strings(2*numel(subarrayLengths),1);
    hIndex = 0;

    for i = 1:numel(subarrayLengths)
        L = subarrayLengths(i);

        direct = T( ...
            T.Method=="SS-MV Direct" & ...
            T.SubarrayLength_L==L,:);
        direct = sortrows(direct,'NumSensors_M');

        sm = T( ...
            T.Method=="SM-MV" & ...
            T.SubarrayLength_L==L,:);
        sm = sortrows(sm,'NumSensors_M');

        hIndex = hIndex+1;
        runtimeHandles(hIndex) = errorbar( ...
            ax1, ...
            direct.NumSensors_M, ...
            1000*direct.MedianRuntime_Sec, ...
            500*direct.RuntimeIQR_Sec, ...
            '-o', ...
            'Color',colors(i,:), ...
            'MarkerFaceColor','w', ...
            'MarkerSize',6, ...
            'LineWidth',1.4, ...
            'CapSize',6);
        runtimeLabels(hIndex) = ...
            "Direct, L="+string(L);

        hIndex = hIndex+1;
        runtimeHandles(hIndex) = errorbar( ...
            ax1, ...
            sm.NumSensors_M, ...
            1000*sm.MedianRuntime_Sec, ...
            500*sm.RuntimeIQR_Sec, ...
            '--s', ...
            'Color',colors(i,:), ...
            'MarkerFaceColor',colors(i,:), ...
            'MarkerSize',6, ...
            'LineWidth',1.4, ...
            'CapSize',6);
        runtimeLabels(hIndex) = ...
            "SM, L="+string(L);
    end

    set(ax1,'XScale','log','YScale','log');
    ax1.XTick = [32 64 128 256];
    ax1.XTickLabel = {'32','64','128','256'};

    xlabel(ax1,'Number of receivers, M', ...
        'FontWeight','bold');
    ylabel(ax1,'Median runtime (ms)', ...
        'FontWeight','bold');
    title(ax1,'(a) Runtime scaling', ...
        'FontWeight','bold');

    legend(ax1,runtimeHandles,runtimeLabels, ...
        'Location','northwest', ...
        'Box','off', ...
        'FontSize',9);

    % ------------------------------------------------------------
    % Panel (b): slopes
    % ------------------------------------------------------------
    ax2 = nexttile(tl,2);
    hold(ax2,'on');
    box(ax2,'on');
    grid(ax2,'on');

    directSlopes = NaN(numel(subarrayLengths),1);
    smSlopes = NaN(numel(subarrayLengths),1);
    directR2 = NaN(numel(subarrayLengths),1);
    smR2 = NaN(numel(subarrayLengths),1);

    for i = 1:numel(subarrayLengths)
        L = subarrayLengths(i);

        rowD = slopes( ...
            slopes.Method=="SS-MV Direct" & ...
            slopes.SubarrayLength_L==L,:);

        rowS = slopes( ...
            slopes.Method=="SM-MV" & ...
            slopes.SubarrayLength_L==L,:);

        if ~isempty(rowD)
            directSlopes(i) = rowD.LogLogSlope(1);
            directR2(i) = rowD.RSquared(1);
        end

        if ~isempty(rowS)
            smSlopes(i) = rowS.LogLogSlope(1);
            smR2(i) = rowS.RSquared(1);
        end
    end

    B = bar(ax2, ...
        categorical(string(subarrayLengths)), ...
        [directSlopes smSlopes], ...
        'grouped');

    B(1).DisplayName = 'Direct SS-MV';
    B(2).DisplayName = 'Sequential SM-MV';

    yline(ax2,1,'--','Linear slope', ...
        'Color',[0.25 0.25 0.25], ...
        'LabelHorizontalAlignment','left');

    ylabel(ax2,'Fitted log-log slope, b', ...
        'FontWeight','bold');
    xlabel(ax2,'Subarray length, L', ...
        'FontWeight','bold');
    title(ax2,'(b) Empirical scaling exponent', ...
        'FontWeight','bold');

    ylim(ax2,[0 max(1.55,1.15*max([directSlopes;smSlopes]))]);

    legend(ax2,'Location','northwest','Box','off');

    for i = 1:numel(subarrayLengths)
        text(ax2,i-0.16,directSlopes(i)+0.045, ...
            sprintf('%.3f\nR^2=%.3f', ...
            directSlopes(i),directR2(i)), ...
            'HorizontalAlignment','center', ...
            'FontSize',8);

        text(ax2,i+0.16,smSlopes(i)+0.045, ...
            sprintf('%.3f\nR^2=%.3f', ...
            smSlopes(i),smR2(i)), ...
            'HorizontalAlignment','center', ...
            'FontSize',8);
    end

    % ------------------------------------------------------------
    % Panel (c): runtime ratio
    % ------------------------------------------------------------
    ax3 = nexttile(tl,3);
    hold(ax3,'on');
    box(ax3,'on');
    grid(ax3,'on');

    ratioHandles = gobjects(numel(subarrayLengths),1);
    ratioLabels = strings(numel(subarrayLengths),1);

    for i = 1:numel(subarrayLengths)
        L = subarrayLengths(i);

        block = R(R.SubarrayLength_L==L,:);
        block = sortrows(block,'NumSensors_M');

        ratioHandles(i) = plot( ...
            ax3, ...
            block.NumSensors_M, ...
            block.SM_to_Direct_Ratio, ...
            '-o', ...
            'Color',colors(i,:), ...
            'MarkerFaceColor',colors(i,:), ...
            'MarkerSize',6, ...
            'LineWidth',1.4);

        ratioLabels(i) = "L="+string(L);
    end

    yline(ax3,1,'--','Equal runtime', ...
        'Color',[0.25 0.25 0.25], ...
        'LabelHorizontalAlignment','left');

    ax3.XScale = 'log';
    ax3.XTick = [32 64 128 256];
    ax3.XTickLabel = {'32','64','128','256'};

    xlabel(ax3,'Number of receivers, M', ...
        'FontWeight','bold');
    ylabel(ax3,'SM-MV / direct runtime', ...
        'FontWeight','bold');
    title(ax3,'(c) Relative runtime', ...
        'FontWeight','bold');

    legend(ax3,ratioHandles,ratioLabels, ...
        'Location','northwest', ...
        'Box','off');

    allAxes = [ax1 ax2 ax3];
    for ax = allAxes
        ax.FontName = 'Times New Roman';
        ax.FontSize = 11;
        ax.LineWidth = 0.9;
        ax.GridAlpha = 0.20;
    end

    title(tl, ...
        'Computational scaling of direct SS-MV and sequential SM-MV', ...
        'FontName','Times New Roman', ...
        'FontSize',18, ...
        'FontWeight','bold');

    files = exportFigureSet( ...
        fig,outputDir,'figure_mv_complexity');

    close(fig);
end

% ========================================================================
function files = exportFigureSet(fig,outputDir,baseName)

    pngFile = fullfile(outputDir,baseName + ".png");
    pdfFile = fullfile(outputDir,baseName + ".pdf");
    figFile = fullfile(outputDir,baseName + ".fig");

    drawnow;

    exportgraphics(fig,pngFile,'Resolution',300);
    exportgraphics(fig,pdfFile,'ContentType','vector');
    savefig(fig,figFile);

    files = [string(pngFile);string(pdfFile);string(figFile)];
end

% ========================================================================
function latexFile = writeLatexSnippet(outputDir)

    latexFile = fullfile(outputDir,'latex_include_snippet.tex');
    fid = fopen(latexFile,'w');

    if fid < 0
        warning('Could not create LaTeX snippet.');
        return;
    end

    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid,'%% Add to preamble:\n');
    fprintf(fid,'%% \\\\usepackage{graphicx}\n\n');

    fprintf(fid,'\\\\begin{figure*}[t]\n');
    fprintf(fid,'\\\\centering\n');
    fprintf(fid,'\\\\includegraphics[width=0.98\\\\textwidth]{figure_resolution.pdf}\n');
    fprintf(fid,'\\\\caption{Point-target resolution and peak-sidelobe behavior at 20 and 0~dB. Symbols show the mean over two target positions and five noise realizations; error bars denote one standard deviation. Full-view coherence and variance weighting generally narrowed the measured main lobe and reduced sidelobe levels relative to DAS. Under limited-view acquisition, radial FWHM was reduced by SM-MV and its weighted variants, whereas tangential sidelobes remained pronounced. FWHM narrowing after nonlinear weighting should therefore be interpreted together with extended-target fidelity and branch preservation.}\n');
    fprintf(fid,'\\\\label{fig:resolution}\n');
    fprintf(fid,'\\\\end{figure*}\n\n');

    fprintf(fid,'\\\\begin{figure*}[t]\n');
    fprintf(fid,'\\\\centering\n');
    fprintf(fid,'\\\\includegraphics[width=0.98\\\\textwidth]{figure_mv_complexity.pdf}\n');
    fprintf(fid,'\\\\caption{Computational scaling of direct SS-MV and sequential SM-MV. Panel (a) shows median runtime over seven timed repetitions; centered error bars have total length equal to the reported runtime IQR because the source table does not contain separate quartiles. Panel (b) shows fitted slopes from $\\\\log T=a+b\\\\log M$, and panel (c) shows the runtime ratio. SM-MV exhibited approximately linear to mildly superlinear scaling for fixed subarray length but remained slower than the optimized direct MATLAB solver over the tested range. Numerical equivalence is reported separately from the dedicated equivalence experiments.}\n');
    fprintf(fid,'\\\\label{fig:mv_complexity}\n');
    fprintf(fid,'\\\\end{figure*}\n');
end

% ========================================================================
function validateResolutionColumns(T)

    required = { ...
        'SensorGeometry','Algorithm','RequestedSNR_dB', ...
        'FWHM_Radial_mm','FWHM_Tangential_mm', ...
        'PSL_Radial_dB','PSL_Tangential_dB', ...
        'PeakLocalizationError_mm'};

    missing = required(~ismember(required,T.Properties.VariableNames));

    if ~isempty(missing)
        error('resolution_results.csv is missing: %s', ...
            strjoin(missing,', '));
    end
end

% ========================================================================
function validateComplexityColumns(T)

    required = { ...
        'Method','NumSensors_M','SubarrayLength_L', ...
        'MedianRuntime_Sec','RuntimeIQR_Sec'};

    missing = required(~ismember(required,T.Properties.VariableNames));

    if ~isempty(missing)
        error('mv_complexity_results.csv is missing: %s', ...
            strjoin(missing,', '));
    end
end

% ========================================================================
function validateSlopeColumns(T)

    required = { ...
        'Method','SubarrayLength_L', ...
        'LogLogSlope','RSquared'};

    missing = required(~ismember(required,T.Properties.VariableNames));

    if ~isempty(missing)
        error('mv_complexity_slopes.csv is missing: %s', ...
            strjoin(missing,', '));
    end
end

% ========================================================================
function m = meanFinite(x)

    x = x(isfinite(x));

    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end

% ========================================================================
function s = stdFinite(x)

    x = x(isfinite(x));

    if numel(x) < 2
        s = NaN;
    else
        s = std(x,0);
    end
end

% ========================================================================
function assertFile(fileName)

    if ~isfile(fileName)
        error('File not found: %s',fileName);
    end
end

% ========================================================================
function ensureFolder(folderName)

    if ~isfolder(folderName)
        [ok,msg] = mkdir(folderName);

        if ~ok
            error('Unable to create folder %s: %s',folderName,msg);
        end
    end
end
