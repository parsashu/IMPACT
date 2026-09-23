function PACT_AdaptiveCohort()
% PACT_AdaptiveCohort
% -------------------------------------------------------------------------
% adaptive-cohort completion only.
%
% PURPOSE
%   Extend the limited-view adaptive cohort from the original 5 phantoms to
%   all 11 phantoms without rerunning k-Wave and without changing any
%   beamformer parameter.
%
% NEW RECONSTRUCTIONS ONLY
%   Missing phantoms:
%       Circle, Ring, Two_Circles, Ellipse, Rectangle, Triangle
%   Geometry:
%       LIMITED_VIEW_ARC only
%   Algorithms:
%       SM-MV, CF-SM-MV, DASDSF-SM-MV, NCVDAS-SM-MV
%   SNR:
%       Inf, 20, 10, 0, -5 dB
%   Repetitions:
%       Inf -> 1
%       finite SNR -> 5
%
% Expected new rows:
%   6 phantoms * (1 + 4*5 cases) * 4 algorithms = 504 rows.
%
% IMPORTANT
%   * No k-Wave forward simulation is called.
%   * Existing v9.3 outputs are read-only.
%   * Clean cached full-ring sensor traces are reused.
%   * Noise seeds reproduce the original main-benchmark seed rule.
%   * Primary MV / weighting parameters are loaded from the original v9.3
%     simulation_configuration.mat and are NOT retuned.
%   * The run is crash-safe and resume-safe.
%
% Companion analysis:
%   analyze_adaptive_cohort.m
% -------------------------------------------------------------------------

close all;
clc;

projectRoot = pwd;
cacheRoot = fullfile(projectRoot,'benchmark_v93_reviewer_revision');
outputDir = fullfile(projectRoot,'ADAPTIVE_11PHANTOM_COMPLETION');

if ~isfolder(cacheRoot)
    error('PACT:MissingCacheRoot', ...
        ['Expected the completed v9.3 folder here:\n%s\n\n' ...
         'Run this file from the same project folder that contains ' ...
         'benchmark_v93_reviewer_revision.'], cacheRoot);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

configFile = fullfile(cacheRoot,'simulation_configuration.mat');
if ~isfile(configFile)
    error('PACT:MissingConfiguration', ...
        'Missing original configuration file: %s',configFile);
end

loadedCfg = load(configFile,'cfg','background','sensorXFull', ...
    'sensorYFull','limitedIndices','phantomNames');
requiredCfgFields = {'cfg','background','sensorXFull','sensorYFull', ...
    'limitedIndices','phantomNames'};
for k = 1:numel(requiredCfgFields)
    if ~isfield(loadedCfg,requiredCfgFields{k})
        error('PACT:IncompleteConfiguration', ...
            'simulation_configuration.mat is missing field "%s".', ...
            requiredCfgFields{k});
    end
end

cfg = loadedCfg.cfg;
background = loadedCfg.background;
sensorXFull = loadedCfg.sensorXFull;
sensorYFull = loadedCfg.sensorYFull;
limitedIndices = loadedCfg.limitedIndices;
phantomNames = string(loadedCfg.phantomNames);

% Keep the original benchmark parameters exactly.
originalOutputDir = cfg.outputDir; %#ok<NASGU>
cfg.outputDir = cacheRoot;
cfg.resume = true;

targetPhantoms = ["Circle","Ring","Two_Circles","Ellipse","Rectangle","Triangle"];
geometry = "LIMITED_VIEW_ARC";
algorithms = ["SM-MV","CF-SM-MV","DASDSF-SM-MV","NCVDAS-SM-MV"];
snrValues = [Inf,20,10,0,-5];
finiteRealizations = 5;

% These checks prevent accidental algorithm retuning in the reviewer run.
fprintf('===============================================================\n');
fprintf(' PACT v9.6 - ADAPTIVE 11-PHANTOM COMPLETION ONLY\n');
fprintf('===============================================================\n');
fprintf('Cache root : %s\n',cacheRoot);
fprintf('Output     : %s\n',outputDir);
fprintf('Geometry   : %s\n',geometry);
fprintf('Algorithms : %s\n',strjoin(cellstr(algorithms),', '));
fprintf('MV L       : %d\n',cfg.mvSubarrayLength);
fprintf('MV loading : %.6g\n',cfg.mvDiagonalLoading);
fprintf('DASDSF q   : %.6g\n',cfg.dasdsfExponent);
fprintf('NCVDAS q   : %.6g\n',cfg.ncVdasExponent);
fprintf('NCVDAS beta: %.6g\n',cfg.ncVdasNoisePenalty);
fprintf('No k-Wave forward simulation will be executed.\n');
fprintf('Expected new result rows: 504\n');
fprintf('===============================================================\n\n');

% Reconstruct the exact benchmark coordinate system and straight-ray delays.
xVec = ((1:cfg.Nx) - (cfg.Nx + 1)/2) * cfg.dx;
yVec = ((1:cfg.Ny) - (cfg.Ny + 1)/2) * cfg.dy;
[Y,X] = meshgrid(yVec,xVec);
tArray = (0:cfg.Nt-1)*cfg.dt;
roiMask = background.bodyMask;

fprintf('Precomputing %s travel-time maps once...\n',char(cfg.delayModel));
travelTimesFull = computeTravelTimes(X,Y,xVec,yVec, ...
    sensorXFull,sensorYFull,background.soundSpeed,cfg);
sensorIndices = limitedIndices;
travelTimes = travelTimesFull(sensorIndices,:);
numSensors = numel(sensorIndices);

% Strict preflight: every cache and evaluation file must already exist.
for phantomName = targetPhantoms
    phantomDir = fullfile(cacheRoot,char(phantomName));
    cleanFile = fullfile(phantomDir,'cached_forward_data','clean_full_ring_data.mat');
    physicalFile = fullfile(phantomDir,'physical_maps.mat');
    maskFile = fullfile(phantomDir,'evaluation_masks.mat');
    if ~isfile(cleanFile)
        error('PACT:MissingForwardCache', ...
            ['Missing cached forward data for %s:\n%s\n' ...
             'This runner will NOT launch k-Wave.'],char(phantomName),cleanFile);
    end
    if ~isfile(physicalFile)
        error('PACT:MissingPhysicalMaps','Missing: %s',physicalFile);
    end
    if ~isfile(maskFile)
        error('PACT:MissingEvaluationMasks','Missing: %s',maskFile);
    end
    if ~any(phantomNames==phantomName)
        error('PACT:UnknownPhantom','%s is absent from the original phantom list.', ...
            char(phantomName));
    end
end
fprintf('Cache/evaluation preflight: PASS for all six missing phantoms.\n');

% Verify that DAS reference rows already exist for the full 11-phantom
% adaptive schedule. They are needed only for the companion effect-size
% analysis; DAS is not reconstructed here.
mainCsv = fullfile(cacheRoot,'all_reconstruction_results_PARTIAL.csv');
if ~isfile(mainCsv)
    mainCsv = fullfile(cacheRoot,'all_reconstruction_results.csv');
end
if ~isfile(mainCsv)
    error('PACT:MissingMainResults', ...
        'Could not find the v9.3 main reconstruction CSV in %s.',cacheRoot);
end
mainResults = readtable(mainCsv,'TextType','string');
allExpectedPhantoms = ["Circle","Ring","Two_Circles","Ellipse","Rectangle", ...
    "Triangle","Shepp_Logan","Complex","Medical","Letter_A","Vascular_Tree"];
for phantomName = allExpectedPhantoms
    for snrDb = snrValues
        reps = finiteRealizations;
        if isinf(snrDb), reps = 1; end
        for realization = 1:reps
            idx = mainResults.Phantom==phantomName & ...
                mainResults.SensorGeometry==geometry & ...
                mainResults.Algorithm=="DAS" & ...
                snrVectorEqual(mainResults.RequestedSNR_dB,snrDb) & ...
                mainResults.NoiseRealization==realization;
            if nnz(idx)~=1
                error('PACT:MissingDASReference', ...
                    ['Expected exactly one DAS reference row for %s, %s, ' ...
                     'SNR=%s, realization=%d; found %d.'], ...
                    char(phantomName),char(geometry),formatSNR(snrDb), ...
                    realization,nnz(idx));
            end
        end
    end
end
fprintf('DAS reference preflight: PASS for all 11 phantoms.\n\n');

partialCsv = fullfile(outputDir,'adaptive_11phantom_results_PARTIAL.csv');
partialMat = fullfile(outputDir,'adaptive_11phantom_results_PARTIAL.mat');
finalCsv = fullfile(outputDir,'adaptive_11phantom_results.csv');
finalMat = fullfile(outputDir,'adaptive_11phantom_results.mat');

if isfile(partialCsv)
    results = readtable(partialCsv,'TextType','string');
    fprintf('Resume mode: loaded %d existing rows.\n',height(results));
else
    results = emptyAdaptiveCompletionTable();
end
results = uniqueAdaptiveCompletionRows(results);

for phantomCounter = 1:numel(targetPhantoms)
    phantomName = targetPhantoms(phantomCounter);
    p = find(phantomNames==phantomName,1);

    phantomDir = fullfile(cacheRoot,char(phantomName));
    loadedPhysical = load(fullfile(phantomDir,'physical_maps.mat'),'optical');
    loadedMasks = load(fullfile(phantomDir,'evaluation_masks.mat'), ...
        'roiMask','targetMask','backgroundMask','gtNorm');

    optical = loadedPhysical.optical;
    gtNorm = loadedMasks.gtNorm;
    localRoiMask = loadedMasks.roiMask;
    targetMask = loadedMasks.targetMask;
    backgroundMask = loadedMasks.backgroundMask;

    loadedClean = load(fullfile(phantomDir,'cached_forward_data', ...
        'clean_full_ring_data.mat'),'cleanFullData');
    cleanData = double(loadedClean.cleanFullData(sensorIndices,:));

    fprintf('\n[%d/%d] %s | limited view | %d sensors\n', ...
        phantomCounter,numel(targetPhantoms),char(phantomName),numSensors);

    for snrDb = snrValues
        reps = finiteRealizations;
        if isinf(snrDb), reps = 1; end

        for realization = 1:reps
            if allAdaptiveCompletionRowsExist(results,phantomName,geometry, ...
                    snrDb,realization,algorithms)
                fprintf('  SNR %7s | rep %d/%d : already complete\n', ...
                    formatSNR(snrDb),realization,reps);
                continue;
            end

            % EXACT same seed rule as the original v9.3 main benchmark.
            rng(caseSeed(cfg.masterSeed,p,geometry,snrDb,realization),'twister');
            [noisyData,knownNoiseSigma,measuredSNR] = ...
                addNoiseAtRequestedSNR(cleanData,snrDb);
            estimatedNoiseVariance = estimatePretriggerNoiseVariance( ...
                noisyData,cfg.noiseEstimationSamples);
            estimatedNoiseSigma = sqrt(max(estimatedNoiseVariance,0));

            tic;
            delayed = sampleDelayedData(noisyData,tArray,travelTimes);
            delayRuntime = toc;

            context.isCircularArray = false;
            context.pixelMask = localRoiMask;
            context.noiseVariance = estimatedNoiseVariance;
            context.numSensors = numSensors;

            outputs = reconstructAlgorithmSet(delayed,algorithms,cfg,context);

            for a = 1:numel(outputs)
                algorithm = string(outputs(a).Algorithm);
                if adaptiveCompletionRowExists(results,phantomName,geometry, ...
                        snrDb,realization,algorithm)
                    continue;
                end

                recNorm = normalizeAmplitude(outputs(a).RawImage,localRoiMask);
                recNorm(~localRoiMask) = 0;
                metrics = computeArticleMetrics(recNorm,gtNorm,localRoiMask, ...
                    targetMask,backgroundMask,cfg.segmentationThreshold, ...
                    cfg.dx,cfg.gCNRBins);

                totalRuntime = delayRuntime + outputs(a).BeamformRuntimeSec;
                row = makeAdaptiveCompletionRow(phantomName,geometry,algorithm, ...
                    snrDb,measuredSNR,realization,knownNoiseSigma, ...
                    estimatedNoiseSigma,metrics,delayRuntime, ...
                    outputs(a).BeamformRuntimeSec,totalRuntime,numSensors,cfg);
                results = [results;row]; %#ok<AGROW>

                fprintf(['    %-16s | SNR %7s | rep %d | PSNR %6.2f | ' ...
                    'SSIM %.4f | Dice %.3f | gCNR %.3f\n'], ...
                    char(algorithm),formatSNR(snrDb),realization, ...
                    metrics.PSNR,metrics.SSIM,metrics.Dice,metrics.gCNR);
            end

            results = uniqueAdaptiveCompletionRows(results);
            writetable(results,partialCsv);
            save(partialMat,'results','cfg','targetPhantoms','algorithms', ...
                'snrValues','finiteRealizations','-v7.3');

            fprintf('  checkpoint -> %d / 504 rows\n',height(results));
        end
    end
end

results = uniqueAdaptiveCompletionRows(results);
expectedRows = numel(targetPhantoms) * (1 + 4*finiteRealizations) * numel(algorithms);
if height(results) ~= expectedRows
    error('PACT:IncompleteAdaptiveCompletion', ...
        'Expected %d rows but obtained %d.',expectedRows,height(results));
end

writetable(results,finalCsv);
save(finalMat,'results','cfg','targetPhantoms','algorithms', ...
    'snrValues','finiteRealizations','-v7.3');

% Small provenance record.
fid = fopen(fullfile(outputDir,'ADAPTIVE_11PHANTOM_RUN_MANIFEST.txt'),'w');
if fid>=0
    fprintf(fid,'PACT v9.6 adaptive-only completion\n');
    fprintf(fid,'Source cache root: %s\n',cacheRoot);
    fprintf(fid,'New phantoms: %s\n',strjoin(cellstr(targetPhantoms),', '));
    fprintf(fid,'Geometry: %s\n',geometry);
    fprintf(fid,'Algorithms: %s\n',strjoin(cellstr(algorithms),', '));
    fprintf(fid,'SNR values: Inf, 20, 10, 0, -5 dB\n');
    fprintf(fid,'Finite-SNR realizations: %d\n',finiteRealizations);
    fprintf(fid,'Expected/final rows: %d\n',height(results));
    fprintf(fid,'Master seed: %d\n',cfg.masterSeed);
    fprintf(fid,'MV subarray length: %d\n',cfg.mvSubarrayLength);
    fprintf(fid,'MV diagonal loading: %.16g\n',cfg.mvDiagonalLoading);
    fprintf(fid,'DASDSF exponent: %.16g\n',cfg.dasdsfExponent);
    fprintf(fid,'NCVDAS exponent: %.16g\n',cfg.ncVdasExponent);
    fprintf(fid,'NCVDAS noise penalty: %.16g\n',cfg.ncVdasNoisePenalty);
    fprintf(fid,'Forward simulations executed by this file: 0\n');
    fclose(fid);
end

fprintf('\n===============================================================\n');
fprintf(' ADAPTIVE COMPLETION FINISHED\n');
fprintf(' Final rows : %d\n',height(results));
fprintf(' CSV        : %s\n',finalCsv);
fprintf(' MAT        : %s\n',finalMat);
fprintf('===============================================================\n');

if exist('analyze_adaptive_cohort','file')==2
    fprintf('\nRunning companion 11-phantom effect-size analysis...\n');
    try
        analyze_adaptive_cohort( ...
            finalCsv,mainCsv,fullfile(projectRoot,'ADAPTIVE_11PHANTOM_ANALYSIS'));
    catch ME
        warning('PACT:AdaptiveAnalysisFailed', ...
            ['The 504 reconstruction rows were saved successfully, but the ' ...
             'post-processing step failed: %s'],ME.message);
    end
end
end

%% ========================================================================
function T = emptyAdaptiveCompletionTable()
names = {'Phantom','SensorGeometry','Algorithm','RequestedSNR_dB', ...
    'MeasuredSNR_dB','NoiseRealization','KnownNoiseSigma_Pa', ...
    'EstimatedNoiseSigma_Pa','RMSE_ROI','PSNR_ROI_dB','SSIM','Dice', ...
    'Hausdorff_mm','HD95_mm','CNR','gCNR','TargetSNR','Contrast_dB', ...
    'ArtifactToTarget_dB','TargetMean','BackgroundMean','TargetStd', ...
    'BackgroundStd','DelayRuntime_Sec','BeamformRuntime_Sec','TotalRuntime_Sec', ...
    'NumSensors','MVSubarrayLength','MVDiagonalLoading','DASDSFExponent', ...
    'NCVDASExponent','NCVDASNoisePenalty','NoiseVarianceSource','RunTag'};
types = [repmat({'string'},1,3),repmat({'double'},1,29),{'string','string'}];
T = table('Size',[0,numel(names)],'VariableTypes',types,'VariableNames',names);
end

%% ========================================================================
function row = makeAdaptiveCompletionRow(phantomName,geometry,algorithm, ...
    requestedSNR,measuredSNR,realization,knownNoiseSigma,estimatedNoiseSigma, ...
    metrics,delayRuntime,beamformRuntime,totalRuntime,numSensors,cfg)
row = emptyAdaptiveCompletionTable();
row(1,:) = {string(phantomName),string(geometry),string(algorithm), ...
    requestedSNR,measuredSNR,realization,knownNoiseSigma,estimatedNoiseSigma, ...
    metrics.RMSE,metrics.PSNR,metrics.SSIM,metrics.Dice,metrics.HausdorffMM, ...
    metrics.HD95MM,metrics.CNR,metrics.gCNR,metrics.TargetSNR, ...
    metrics.ContrastDB,metrics.ArtifactToTargetDB,metrics.TargetMean, ...
    metrics.BackgroundMean,metrics.TargetStd,metrics.BackgroundStd, ...
    delayRuntime,beamformRuntime,totalRuntime,numSensors,cfg.mvSubarrayLength, ...
    cfg.mvDiagonalLoading,cfg.dasdsfExponent,cfg.ncVdasExponent, ...
    cfg.ncVdasNoisePenalty,"PretriggerMAD","V96_ADAPTIVE_11PHANTOM"};
end

%% ========================================================================
function tf = adaptiveCompletionRowExists(T,phantom,geometry,snrDb,realization,algorithm)
if isempty(T), tf = false; return; end
tf = any(T.Phantom==string(phantom) & T.SensorGeometry==string(geometry) & ...
    snrVectorEqual(T.RequestedSNR_dB,snrDb) & ...
    T.NoiseRealization==realization & T.Algorithm==string(algorithm));
end

%% ========================================================================
function tf = allAdaptiveCompletionRowsExist(T,phantom,geometry,snrDb,realization,algorithms)
tf = true;
for a = 1:numel(algorithms)
    tf = tf && adaptiveCompletionRowExists(T,phantom,geometry,snrDb, ...
        realization,algorithms(a));
end
end

%% ========================================================================
function T = uniqueAdaptiveCompletionRows(T)
if isempty(T), return; end
keys = T.Phantom+"|"+T.SensorGeometry+"|"+T.Algorithm+"|"+ ...
    string(T.RequestedSNR_dB)+"|"+string(T.NoiseRealization);
[~,idx] = unique(keys,'stable');
T = T(idx,:);
end

%% ========================================================================
function cfg = defaultConfiguration()
cfg.profile = "PAPER_FULL"; % SMOKE_TEST | PAPER_CORE | PAPER_FULL
cfg.masterSeed = 20260711;
cfg.resume = true;
% Optional: set this to the previous PAPER_FULL output directory to reuse
% clean heterogeneous k-Wave traces and avoid rerunning the forward solver.
% Example: cfg.reuseForwardCacheRoot = 'benchmark_v90_paper_ready';
cfg.reuseForwardCacheRoot = "";

% Grid
cfg.Nx = 256;
cfg.Ny = 256;
cfg.dx = 100e-6;
cfg.dy = 100e-6;
cfg.dt = 18e-9;
cfg.Nt = 1024;

% Geometry
cfg.bodyRadius = 9.0e-3;
cfg.sensorRadius = 10.5e-3;
cfg.numSensorsFull = 256;
cfg.numSensorsLimited = 128;
cfg.elementWidth = 0.20e-3;
cfg.arrayBLITolerance = 0.05;
cfg.arrayUpsamplingRate = 8;
cfg.sensorGeometries = ["CIRCULAR_FULL", "LIMITED_VIEW_ARC"];

% Detector/electronics
cfg.sensorCenterFrequency = 2.0e6;
cfg.sensorBandwidthPercent = 80;

% Acoustic absorption
cfg.alphaPower = 1.10;

% Optical model
cfg.surfaceFluenceJm2 = 100;
cfg.numIlluminationFibers = 8;
cfg.illuminationAngularSigmaDeg = 10;
cfg.illuminationFloorFraction = 0.08;
cfg.targetMuAContrast = 85;
cfg.targetMuSpFraction = 0.12;
cfg.targetGruneisenFraction = 0.08;
cfg.fluenceModel = "DIFFUSION";

% Reconstruction delay model
cfg.delayModel = "STRAIGHT_RAY_HETEROGENEOUS";
cfg.straightRaySamples = 16;
cfg.homogeneousReconSoundSpeed = 1540;

% Noise and repeated experiments
cfg.snrDb = [Inf, 30, 20, 10, 5, 0, -5];
cfg.numNoiseRealizations = 20;
cfg.noiseEstimationSamples = 32;

% Algorithms
cfg.cheapAlgorithms = ["DAS", "CF-DAS", "SCF-DAS", "DMAS", ...
    "DASDSF", "NC-CF", "NC-VDAS"];
cfg.adaptiveAlgorithms = ["SM-MV", "CF-SM-MV", ...
    "DASDSF-SM-MV", "NCVDAS-SM-MV"];
cfg.mainPhantomNames = "ALL";
cfg.adaptivePhantomNames = ["Shepp_Logan", "Complex", "Medical", ...
    "Letter_A", "Vascular_Tree"];
cfg.adaptiveSNRdB = [Inf, 20, 10, 0, -5];
cfg.adaptiveNoiseRealizations = 5;

% Coherence/variance parameters
cfg.scfExponent = 1.5;
cfg.dasdsfExponent = 1.0;
cfg.ncVdasExponent = 1.0;
cfg.ncVdasNoisePenalty = 1.0;

% MV parameters
cfg.mvSubarrayLength = 16;
cfg.mvDiagonalLoading = 1e-2;
cfg.mvLoadingFloorRelative = 1e-10;
cfg.mvTilePixels = 64;
cfg.mvValidationPhantomNames = ["Medical", "Vascular_Tree"];
cfg.mvValidationSNRdB = [Inf, 10, 0];
cfg.mvValidationRealization = 1;

% Evaluation
cfg.segmentationThreshold = 0.30;
cfg.backgroundGtMaximum = 0.05;
cfg.backgroundMarginMM = 0.60;
cfg.gCNRBins = 128;

% Saved reconstructions
cfg.saveSNRdB = [20, 0, -5];
cfg.saveRealization = 1;
cfg.saveRawReconstructions = true;

% Resolution benchmark
cfg.runResolutionBenchmark = true;
cfg.resolutionSNRdB = [20, 0];
cfg.resolutionRealizations = 5;
cfg.resolutionPointSigmaMM = 0.12;
cfg.resolutionPointLocationsMM = [ ...
    0.0,  0.0; ...
    3.0,  0.0;  0.0,  3.0; -3.0,  0.0;  0.0, -3.0; ...
    4.24, 4.24; -4.24, 4.24; -4.24,-4.24; 4.24,-4.24];
cfg.resolutionAlgorithms = [cfg.cheapAlgorithms, cfg.adaptiveAlgorithms];
cfg.profileHalfWidthMM = 4.0;
cfg.profileSamples = 1601;

% Complexity/equivalence benchmark
cfg.runComplexityBenchmark = true;
cfg.complexitySensorCounts = [32, 64, 128, 256];
cfg.complexitySubarrayLengths = [8, 16, 32];
cfg.complexityNumPixels = 256;
cfg.complexityRepetitions = 7;
cfg.complexityWarmups = 1;
cfg.inverseResidualSamples = 8;

% Reviewer-driven parameter sensitivity (predeclared; no outcome-based tuning)
cfg.runParameterSensitivity = true;
cfg.sensitivityPhantomNames = ["Medical","Vascular_Tree","Complex"];
cfg.sensitivityGeometries = ["CIRCULAR_FULL","LIMITED_VIEW_ARC"];
cfg.sensitivitySNRdB = [20,0,-5];
cfg.sensitivityRealizations = 1:3;
cfg.sensitivitySCFExponent = [0.5, 1.0, 1.5, 2.0, 3.0];
cfg.sensitivityDASDSFExponent = [0.5, 1.0, 1.5, 2.0];
cfg.sensitivityNCVDASExponent = [0.5, 1.0, 1.5, 2.0];
cfg.sensitivityNCVDASPenalty = [0.25, 0.5, 1.0, 2.0];
cfg.sensitivityMVSubarray = [8, 16, 24];
cfg.sensitivityMVLoading = [1e-3, 1e-2, 1e-1];

% Reviewer-driven acoustic-model mismatch benchmark. These conditions reuse
% cached heterogeneous k-Wave data and change ONLY reconstruction delays.
cfg.runDelayMismatchBenchmark = true;
cfg.mismatchPhantomNames = ["Medical","Vascular_Tree","Complex"];
cfg.mismatchGeometries = ["CIRCULAR_FULL","LIMITED_VIEW_ARC"];
cfg.mismatchSNRdB = [20,0];
cfg.mismatchRealizations = 3;
cfg.mismatchAlgorithms = ["DAS","DMAS","DASDSF","NC-CF","NC-VDAS", ...
    "SM-MV","CF-SM-MV","DASDSF-SM-MV","NCVDAS-SM-MV"];
cfg.mismatchConditions = ["BASELINE_HETERO","HOMOGENEOUS_1540", ...
    "GLOBAL_MINUS_2PCT","GLOBAL_PLUS_2PCT","MISREG_X_0P5MM", ...
    "LOWRES_X4","UNDERSEGMENT_2PX"];
cfg.mismatchGlobalFraction = 0.02;
cfg.mismatchShiftMM = 0.50;
cfg.mismatchDownsampleFactor = 4;
cfg.mismatchSegmentationPixels = 2;

% Reviewer-driven noise robustness benchmark. These conditions also reuse
% cached clean sensor data; no new k-Wave forward solve is required.
cfg.runNoiseRobustnessBenchmark = true;
cfg.noiseRobustnessPhantomNames = ["Medical","Vascular_Tree","Complex"];
cfg.noiseRobustnessGeometries = ["CIRCULAR_FULL","LIMITED_VIEW_ARC"];
cfg.noiseRobustnessSNRdB = -5;
cfg.noiseRobustnessRealizations = 5;
cfg.noiseRobustnessAlgorithms = ["DAS","DMAS","DASDSF","NC-CF","NC-VDAS", ...
    "SM-MV","CF-SM-MV","DASDSF-SM-MV","NCVDAS-SM-MV"];
cfg.noiseRobustnessModels = ["IID_GAUSSIAN","CHANNEL_CORRELATED_GAUSSIAN", ...
    "HETEROSCEDASTIC_GAUSSIAN","LAPLACE_IID"];
cfg.channelNoiseCorrelation = 0.30;
cfg.channelNoiseScaleSpread = 0.20;

% Statistical analysis
cfg.permutationCount = 10000;
cfg.bootstrapCount = 4000;
cfg.statisticsSeed = 661177;

cfg.outputDir = 'benchmark_v93_reviewer_revision';
end

%% ========================================================================
function cfg = applyRunProfile(cfg, profile)
profile = upper(string(profile));
cfg.profile = profile;

switch profile
    case "SMOKE_TEST"
        cfg.Nx = 128;
        cfg.Ny = 128;
        cfg.dx = 150e-6;
        cfg.dy = cfg.dx;
        cfg.dt = 25e-9;
        cfg.Nt = 768;
        cfg.bodyRadius = 6.0e-3;
        cfg.sensorRadius = 7.5e-3;
        cfg.numSensorsFull = 96;
        cfg.numSensorsLimited = 48;
        cfg.elementWidth = 0.25e-3;
        cfg.snrDb = [Inf, 10];
        cfg.numNoiseRealizations = 1;
        cfg.mainPhantomNames = ["Medical"];
        cfg.adaptivePhantomNames = ["Medical"];
        cfg.adaptiveSNRdB = [Inf, 10];
        cfg.adaptiveNoiseRealizations = 1;
        cfg.mvSubarrayLength = 8;
        cfg.mvTilePixels = 32;
        cfg.runResolutionBenchmark = false;
        cfg.runComplexityBenchmark = true;
        cfg.complexitySensorCounts = [32, 64];
        cfg.complexitySubarrayLengths = [8];
        cfg.complexityNumPixels = 32;
        cfg.complexityRepetitions = 2;
        cfg.runParameterSensitivity = false;
        cfg.runDelayMismatchBenchmark = false;
        cfg.runNoiseRobustnessBenchmark = false;
        cfg.permutationCount = 1000;
        cfg.bootstrapCount = 500;
        cfg.outputDir = 'benchmark_v93_smoke_test';

    case "PAPER_CORE"
        cfg.snrDb = [Inf, 20, 10, 0];
        cfg.numNoiseRealizations = 5;
        cfg.mainPhantomNames = ["Circle", "Shepp_Logan", "Complex", ...
            "Medical", "Vascular_Tree"];
        cfg.adaptivePhantomNames = ["Shepp_Logan", "Complex", ...
            "Medical", "Vascular_Tree"];
        cfg.adaptiveSNRdB = [Inf, 10, 0];
        cfg.adaptiveNoiseRealizations = 3;
        cfg.resolutionRealizations = 3;
        cfg.mismatchRealizations = 1;
        cfg.noiseRobustnessRealizations = 2;
        cfg.complexityNumPixels = 128;
        cfg.complexityRepetitions = 5;
        cfg.permutationCount = 5000;
        cfg.bootstrapCount = 2000;
        cfg.outputDir = 'benchmark_v93_paper_core';

    case "PAPER_FULL"
        % Use the defaults declared above.

    otherwise
        error('Unknown cfg.profile: %s', char(profile));
end
end

%% ========================================================================
function assertRequiredFunctionsV90()
required = {'kWaveGrid', 'kspaceFirstOrder2D', 'kWaveArray', ...
    'ssim', 'bwperim', 'bwdist', 'histcounts', 'pagemtimes', ...
    'imresize', 'imerode', 'strel', 'exportgraphics', 'jsonencode'};
missing = strings(0,1);
for k = 1:numel(required)
    if exist(required{k}, 'file') == 0 && exist(required{k}, 'class') == 0 && ...
            exist(required{k}, 'builtin') == 0
        missing(end+1,1) = string(required{k}); %#ok<AGROW>
    end
end
if ~isempty(missing)
    error('Missing required MATLAB/k-Wave functions: %s', ...
        strjoin(missing, ', '));
end
end

%% ========================================================================
function validateConfigurationV90(cfg)
if cfg.Nx ~= cfg.Ny || abs(cfg.dx-cfg.dy) > eps
    error('The current implementation requires a square isotropic grid.');
end
if cfg.numSensorsLimited > cfg.numSensorsFull/2
    error('Limited-view sensor count cannot exceed half the full ring.');
end
if cfg.elementWidth >= 2*pi*cfg.sensorRadius/cfg.numSensorsFull
    error('Element width must be smaller than the circumferential pitch.');
end
if cfg.bodyRadius >= cfg.sensorRadius
    error('bodyRadius must be smaller than sensorRadius.');
end
if cfg.mvSubarrayLength < 2 || cfg.mvSubarrayLength > cfg.numSensorsLimited
    error('Invalid mvSubarrayLength.');
end
if cfg.mvDiagonalLoading <= 0
    error('mvDiagonalLoading must be positive for Sherman-Morrison updates.');
end
if cfg.mvTilePixels < 1
    error('mvTilePixels must be positive.');
end
if cfg.noiseEstimationSamples < 4 || cfg.noiseEstimationSamples >= cfg.Nt
    error('noiseEstimationSamples must lie between 4 and Nt-1.');
end
if any(cfg.snrDb(~isinf(cfg.snrDb)) > 100)
    error('Unreasonably high finite SNR requested.');
end
maxExpectedSoundSpeed = 1610;
cfl = maxExpectedSoundSpeed*cfg.dt/min(cfg.dx,cfg.dy);
if cfl > 0.30
    error('Expected CFL %.3f exceeds 0.30.', cfl);
end
minExpectedSoundSpeed = 1430;
requiredTime = (cfg.sensorRadius + cfg.bodyRadius)/minExpectedSoundSpeed;
availableTime = (cfg.Nt-1)*cfg.dt;
if availableTime < requiredTime
    error('Time window %.3f us is shorter than required %.3f us.', ...
        1e6*availableTime, 1e6*requiredTime);
end
end

%% ========================================================================
function tf = isAbsolutePathV91(pathName)
pathName = char(pathName);
if ispc
    tf = ~isempty(regexp(pathName, '^[A-Za-z]:[\\/]', 'once')) || ...
        startsWith(pathName, '\\');
else
    tf = startsWith(pathName, filesep);
end
end

%% ========================================================================
function saveCoreResultsCheckpointV91(outputDir, allResults, ...
    allEquivalence, allSensitivity, cfg)
ensureDirectory(outputDir);
allResults = uniqueMainResults(allResults);
allEquivalence = uniqueEquivalenceResults(allEquivalence);
allSensitivity = uniqueSensitivityResults(allSensitivity);
save(fullfile(outputDir, 'MASTER_CHECKPOINT.mat'), ...
    'allResults', 'allEquivalence', 'allSensitivity', 'cfg', '-v7.3');
writetable(allResults, ...
    fullfile(outputDir, 'all_reconstruction_results_PARTIAL.csv'));
writetable(allEquivalence, ...
    fullfile(outputDir, 'mv_image_equivalence_PARTIAL.csv'));
writetable(allSensitivity, ...
    fullfile(outputDir, 'parameter_sensitivity_PARTIAL.csv'));
end

%% ========================================================================
function publishResultsToBaseWorkspaceV91(cfg, allResults, allEquivalence, ...
    allSensitivity, resolutionResults, complexityResults, complexitySlopes)
assignin('base', 'PACT_cfg', cfg);
assignin('base', 'PACT_allResults', allResults);
assignin('base', 'PACT_allEquivalence', allEquivalence);
assignin('base', 'PACT_allSensitivity', allSensitivity);
assignin('base', 'PACT_resolutionResults', resolutionResults);
assignin('base', 'PACT_complexityResults', complexityResults);
assignin('base', 'PACT_complexitySlopes', complexitySlopes);
end

%% ========================================================================
function saveBenchmarkFailureV91(outputDir, benchmarkName, ME)
ensureDirectory(outputDir);
failureDir = fullfile(outputDir, 'failure_logs');
ensureDirectory(failureDir);
timeStamp = datestr(now, 'yyyymmdd_HHMMSS');
fileName = fullfile(failureDir, ...
    sprintf('%s_failure_%s.mat', benchmarkName, timeStamp));
reportText = getReport(ME, 'extended', 'hyperlinks', 'off');
errorIdentifier = ME.identifier;
errorMessage = ME.message;
errorStack = ME.stack;
save(fileName, 'reportText', 'errorIdentifier', 'errorMessage', 'errorStack');
textFile = strrep(fileName, '.mat', '.txt');
fid = fopen(textFile, 'w');
if fid >= 0
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', reportText);
end
end

%% ========================================================================
function closeFigureIfValidV91(fig)
if ~isempty(fig) && isgraphics(fig)
    close(fig);
end
end

%% ========================================================================
function prepareOutputDirectories(cfg)
ensureDirectory(cfg.outputDir);
ensureDirectory(fullfile(cfg.outputDir, 'paper_figures'));
ensureDirectory(fullfile(cfg.outputDir, 'paper_tables'));
ensureDirectory(fullfile(cfg.outputDir, 'resolution_benchmark'));
ensureDirectory(fullfile(cfg.outputDir, 'complexity_benchmark'));
end

%% ========================================================================
function ensureDirectory(pathName)
if ~exist(pathName, 'dir')
    mkdir(pathName);
end
end

%% ========================================================================
function writeRunManifest(cfg)
fileName = fullfile(cfg.outputDir, 'RUN_MANIFEST.txt');
fid = fopen(fileName, 'w');
if fid < 0
    error('Cannot create %s.', fileName);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'PACT Benchmark v9.3 reviewer revision\n');
fprintf(fid, 'Generated: %s\n', char(datetime('now')));
fprintf(fid, 'Profile: %s\n', char(cfg.profile));
fprintf(fid, 'MATLAB version: %s\n', version);
fprintf(fid, 'Computer architecture: %s\n', computer('arch'));
fprintf(fid, 'No iterative reconstruction/convergence experiment is included.\n');
fprintf(fid, 'Primary proposed components: NC-VDAS and sequential SM-MV.\n');
fprintf(fid, 'Reference MV: direct spatially smoothed, diagonally loaded SS-MV.\n');
fprintf(fid, 'Random seed: %d\n', cfg.masterSeed);
fprintf(fid, 'Output directory: %s\n', cfg.outputDir);
end

%% ========================================================================
function writeConfigurationJSON(cfg)
fileName = fullfile(cfg.outputDir, 'configuration.json');
text = jsonencode(cfg);
fid = fopen(fileName, 'w');
if fid < 0
    warning('Could not write configuration JSON.');
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, text, 'char');
end

%% ========================================================================
function indices = selectPhantomIndices(names, requested)
if isscalar(requested) && upper(string(requested)) == "ALL"
    indices = 1:numel(names);
    return;
end
requested = string(requested);
indices = zeros(1, numel(requested));
for k = 1:numel(requested)
    match = find(names == requested(k), 1);
    if isempty(match)
        error('Requested phantom %s is not available.', char(requested(k)));
    end
    indices(k) = match;
end
end

%% ========================================================================
function [sensorIndices, isCircular] = geometryIndices(geometry, cfg, limitedIndices)
switch upper(string(geometry))
    case "CIRCULAR_FULL"
        sensorIndices = 1:cfg.numSensorsFull;
        isCircular = true;
    case "LIMITED_VIEW_ARC"
        sensorIndices = limitedIndices;
        isCircular = false;
    otherwise
        error('Unknown sensor geometry: %s', char(geometry));
end
end

%% ========================================================================
function algorithms = algorithmsForCase(phantomName, snrDb, realization, cfg)
algorithms = cfg.cheapAlgorithms;
adaptiveSchedule = any(string(phantomName)==cfg.adaptivePhantomNames) && ...
    anySNRMatch(snrDb,cfg.adaptiveSNRdB) && ...
    (isinf(snrDb) || realization<=cfg.adaptiveNoiseRealizations);
if adaptiveSchedule
    algorithms = [algorithms,cfg.adaptiveAlgorithms];
end
algorithms = unique(algorithms,'stable');
end

%% ========================================================================
function tf = isDirectMVValidationCase(phantomName, geometry, snrDb, realization, cfg)
tf = any(string(phantomName) == cfg.mvValidationPhantomNames) && ...
    anySNRMatch(snrDb, cfg.mvValidationSNRdB) && ...
    realization == cfg.mvValidationRealization && ...
    any(string(geometry) == cfg.sensorGeometries);
end

%% ========================================================================
function tf = shouldRunSensitivityCase(phantomName, geometry, snrDb, realization, cfg)
tf = cfg.runParameterSensitivity && ...
    any(string(phantomName) == cfg.sensitivityPhantomNames) && ...
    any(string(geometry) == cfg.sensitivityGeometries) && ...
    anySNRMatch(snrDb, cfg.sensitivitySNRdB) && ...
    any(realization == cfg.sensitivityRealizations);
end

%% ========================================================================
function n = realizationCountForSNR(snrDb, cfg)
if isinf(snrDb)
    n = 1;
else
    n = cfg.numNoiseRealizations;
end
end

%% ========================================================================
function seed = caseSeed(masterSeed, phantomIndex, geometry, snrDb, realization)
geometryCode = sum(double(char(string(geometry))));
if isinf(snrDb)
    snrCode = 999;
else
    snrCode = round(10*(snrDb + 100));
end
seed = mod(masterSeed + 1000003*phantomIndex + 1009*geometryCode + ...
    37*snrCode + realization, 2^31-1);
if seed <= 0
    seed = seed + 12345;
end
end

%% ========================================================================
function [noisyData, sigma, measuredSNR] = addNoiseAtRequestedSNR(cleanData, requestedSNR)
signalRMS = rootMeanSquare(cleanData(:));
if isinf(requestedSNR)
    sigma = 0;
    noisyData = cleanData;
    measuredSNR = Inf;
    return;
end
sigma = signalRMS * 10^(-requestedSNR/20);
noise = sigma * randn(size(cleanData));
noisyData = cleanData + noise;
noiseRMS = rootMeanSquare(noise(:));
measuredSNR = 20*log10(signalRMS/max(noiseRMS, eps));
end

%% ========================================================================
function varianceEstimate = estimatePretriggerNoiseVariance(sensorData, numSamples)
numSamples = min(numSamples, size(sensorData,2));
segment = double(sensorData(:,1:numSamples));
channelMedian = median(segment, 2);
absoluteDeviation = abs(segment - channelMedian);
robustSigma = median(absoluteDeviation, 2) / 0.674489750196082;
classicalSigma = std(segment, 0, 2);
invalid = ~isfinite(robustSigma) | robustSigma <= eps;
robustSigma(invalid) = classicalSigma(invalid);
varianceEstimate = median(robustSigma.^2, 'omitnan');
if isempty(varianceEstimate) || ~isfinite(varianceEstimate) || varianceEstimate < 0
    varianceEstimate = 0;
end
end

%% ========================================================================
function value = rootMeanSquare(x)
x = double(x(:));
value = sqrt(mean(abs(x).^2));
end

%% ========================================================================
function outputs = reconstructAlgorithmSet(delayed, algorithms, cfg, context)
outputs = struct('Algorithm', {}, 'RawImage', {}, 'BeamformRuntimeSec', {});
algorithms = string(algorithms);
needsSM = any(contains(algorithms, "SM-MV"));
needsDirect = any(algorithms == "SS-MV");

smVector = [];
smRuntime = NaN;
if needsSM
    tic;
    smVector = spatiallySmoothedMVTiled(double(delayed), cfg, context, ...
        "SHERMAN_MORRISON");
    smRuntime = toc;
end

directVector = [];
directRuntime = NaN;
if needsDirect
    tic;
    directVector = spatiallySmoothedMVTiled(double(delayed), cfg, context, ...
        "DIRECT");
    directRuntime = toc;
end

for a = 1:numel(algorithms)
    algorithm = algorithms(a);
    switch algorithm
        case "SM-MV"
            vector = smVector;
            runtime = smRuntime;
        case "SS-MV"
            vector = directVector;
            runtime = directRuntime;
        case "CF-SM-MV"
            tic;
            weight = coherenceWeight(double(delayed));
            vector = weight .* smVector;
            runtime = smRuntime + toc;
        case "DASDSF-SM-MV"
            tic;
            weight = dasdsfWeight(double(delayed), cfg.dasdsfExponent);
            vector = weight .* smVector;
            runtime = smRuntime + toc;
        case "NCVDAS-SM-MV"
            tic;
            weight = ncVdasWeight(double(delayed), context.noiseVariance, ...
                cfg.ncVdasExponent, cfg.ncVdasNoisePenalty);
            vector = weight .* smVector;
            runtime = smRuntime + toc;
        otherwise
            tic;
            vector = cheapBeamformerVector(double(delayed), algorithm, ...
                context.noiseVariance, cfg);
            runtime = toc;
    end

    outputs(end+1).Algorithm = algorithm; %#ok<AGROW>
    outputs(end).RawImage = reshape(vector, cfg.Nx, cfg.Ny);
    outputs(end).BeamformRuntimeSec = runtime;
end
end

%% ========================================================================
function y = cheapBeamformerVector(x, algorithm, noiseVariance, cfg)
[M, ~] = size(x);
meanSignal = mean(x, 1);

switch upper(string(algorithm))
    case "DAS"
        y = meanSignal;

    case "CF-DAS"
        w = coherenceWeight(x);
        y = w .* meanSignal;

    case "SCF-DAS"
        signs = sign(x);
        meanSign = mean(signs, 1);
        varianceSign = mean(signs.^2, 1) - meanSign.^2;
        varianceSign = min(max(real(varianceSign), 0), 1);
        w = max(0, 1-sqrt(varianceSign)).^cfg.scfExponent;
        y = w .* meanSignal;

    case "DMAS"
        q = sign(x).*sqrt(abs(x));
        pairSum = 0.5*(sum(q,1).^2 - sum(abs(x),1));
        numPairs = M*(M-1)/2;
        y = pairSum/max(numPairs,1);

    case "DASDSF"
        w = dasdsfWeight(x, cfg.dasdsfExponent);
        y = w .* meanSignal;

    case "NC-CF"
        channelSum = sum(x,1);
        energy = sum(abs(x).^2,1);
        correctedNumerator = max(abs(channelSum).^2 - M*noiseVariance, 0);
        correctedEnergy = max(energy - M*noiseVariance, 0);
        w = correctedNumerator ./ (M*correctedEnergy + eps);
        w = min(max(real(w),0),1);
        y = w .* meanSignal;

    case "NC-VDAS"
        w = ncVdasWeight(x, noiseVariance, cfg.ncVdasExponent, ...
            cfg.ncVdasNoisePenalty);
        y = w .* meanSignal;

    otherwise
        error('Unknown cheap beamformer: %s', char(algorithm));
end
end

%% ========================================================================
function w = coherenceWeight(x)
M = size(x,1);
channelSum = sum(x,1);
energy = sum(abs(x).^2,1);
w = abs(channelSum).^2 ./ (M*energy + eps);
w = min(max(real(w),0),1);
end

%% ========================================================================
function w = dasdsfWeight(x, exponent)
meanSignal = mean(x,1);
centered = x - meanSignal;
stdPopulation = sqrt(mean(abs(centered).^2,1));
ratio = abs(meanSignal)./(stdPopulation + eps);
% Bounded regularization of the DAS-to-standard-deviation ratio.
w = (ratio./(1+ratio)).^exponent;
w = min(max(real(w),0),1);
end

%% ========================================================================
function w = ncVdasWeight(x, noiseVariance, exponent, noisePenalty)
M = size(x,1);
meanSignal = mean(x,1);
populationVariance = mean(abs(x-meanSignal).^2,1);
coherentPower = max(abs(meanSignal).^2 - noiseVariance/M, 0);
incoherentPower = max(populationVariance - (1-1/M)*noiseVariance, 0);
denominator = coherentPower + incoherentPower + ...
    noisePenalty*noiseVariance + eps;
w = (coherentPower./denominator).^exponent;
w = min(max(real(w),0),1);
end

%% ========================================================================
function y = spatiallySmoothedMVTiled(x, cfg, context, method)
[M, numPixels] = size(x);
L = min(cfg.mvSubarrayLength, floor(M/2));
if L < 2
    error('MV subarray length must be at least 2.');
end

if context.isCircularArray
    xExtended = [x; x(1:L-1,:)];
    K = M;
else
    xExtended = x;
    K = M-L+1;
end

activePixels = find(context.pixelMask(:)).';
y = zeros(1, numPixels);
if isempty(activePixels)
    return;
end

globalPower = mean(abs(x(:)).^2);
tileSize = cfg.mvTilePixels;
for first = 1:tileSize:numel(activePixels)
    tilePixels = activePixels(first:min(first+tileSize-1, numel(activePixels)));
    Z = buildSubarrayTensor(xExtended, tilePixels, L, K);
    y(tilePixels) = solveMVTile(Z, context.noiseVariance, ...
        cfg.mvDiagonalLoading, cfg.mvLoadingFloorRelative, ...
        globalPower, method);
end
end

%% ========================================================================
function Z = buildSubarrayTensor(xExtended, pixelIndices, L, K)
B = numel(pixelIndices);
Z = zeros(L, K, B, 'like', xExtended);
for k = 1:K
    Z(:,k,:) = reshape(xExtended(k:k+L-1, pixelIndices), L, 1, B);
end
end

%% ========================================================================
function y = solveMVTile(Z, noiseVariance, loadingFactor, floorRelative, ...
    globalPower, method)
[L,K,B] = size(Z);
meanZ = reshape(mean(Z,2), L, B);
traceR = reshape(sum(sum(abs(Z).^2,1),2), 1, B)/K;
loading = loadingFactor*(traceR/L + noiseVariance) + ...
    floorRelative*max(globalPower, eps);
loading = max(real(loading), eps);
a = ones(L,1);
y = zeros(1,B);

switch upper(string(method))
    case "DIRECT"
        identity = eye(L);
        for b = 1:B
            Zb = Z(:,:,b);
            R = (Zb*Zb')/K + loading(b)*identity;
            u = R\a;
            denominator = a'*u;
            if ~isfinite(denominator) || abs(denominator) <= eps
                y(b) = 0;
            else
                w = u/denominator;
                y(b) = w'*meanZ(:,b);
            end
        end

    case "SHERMAN_MORRISON"
        identityPages = repmat(eye(L), 1, 1, B);
        P = identityPages ./ reshape(loading,1,1,B);
        rootK = sqrt(K);
        for k = 1:K
            U = Z(:,k,:)/rootK;
            PU = pagemtimes(P,U);
            denominator = 1 + sum(conj(U).*PU,1);
            denominator = max(real(denominator), eps);
            correction = pagemtimes(PU, permute(conj(PU), [2,1,3]));
            P = P - correction./denominator;
        end
        aPages = repmat(a,1,1,B);
        Pa = pagemtimes(P,aPages);
        weightDenominator = sum(conj(aPages).*Pa,1);
        weightDenominator = max(real(weightDenominator),eps);
        W = Pa./weightDenominator;
        meanZPages = reshape(meanZ,L,1,B);
        y = reshape(sum(conj(W).*meanZPages,1),1,B);

    otherwise
        error('Unknown MV solver: %s', char(method));
end

y = real(y);
end

%% ========================================================================
function [targetMask, backgroundMask] = makeEvaluationMasks(gt, roiMask, ...
    threshold, backgroundGtMaximum, backgroundMarginMM, pixelSize)
targetMask = roiMask & gt >= threshold;
if ~any(targetMask(:))
    error('Target mask is empty at threshold %.3f.', threshold);
end
marginPixels = backgroundMarginMM/(pixelSize*1e3);
distanceFromTarget = bwdist(targetMask);
backgroundMask = false(size(gt));
thresholdCandidates = unique([backgroundGtMaximum,0.08,0.10,0.15,0.20],'stable');
marginCandidates = unique([marginPixels,max(1,marginPixels/2),1],'stable');
for margin = marginCandidates
    for gtMaximum = thresholdCandidates
        candidate = roiMask & ~targetMask & gt<=gtMaximum & ...
            distanceFromTarget>=margin;
        if nnz(candidate)>=100
            backgroundMask = candidate;
            return;
        end
        if nnz(candidate)>nnz(backgroundMask)
            backgroundMask = candidate;
        end
    end
end
if nnz(backgroundMask)<20
    error('Unable to define a sufficiently large background ROI.');
end
warning('Background ROI contains only %d pixels.',nnz(backgroundMask));
end

%% ========================================================================
function metrics = computeArticleMetrics(rec, gt, roiMask, targetMask, ...
    backgroundMask, threshold, pixelSize, numGCNRBins)
errorValues = rec(roiMask)-gt(roiMask);
metrics.RMSE = sqrt(mean(errorValues.^2));
metrics.PSNR = 20*log10(1/max(metrics.RMSE,eps));

[row,col] = find(roiMask);
rowRange = min(row):max(row);
colRange = min(col):max(col);
metrics.SSIM = ssim(rec(rowRange,colRange), gt(rowRange,colRange), ...
    'DynamicRange',1);

recMask = roiMask & rec >= threshold;
gtMask = roiMask & gt >= threshold;
intersection = nnz(recMask & gtMask);
denominator = nnz(recMask)+nnz(gtMask);
if denominator == 0
    metrics.Dice = 1;
else
    metrics.Dice = 2*intersection/denominator;
end
[metrics.HausdorffMM,metrics.HD95MM] = boundaryDistancesMM( ...
    recMask,gtMask,pixelSize);

targetValues = rec(targetMask);
backgroundValues = rec(backgroundMask);
metrics.TargetMean = mean(targetValues);
metrics.BackgroundMean = mean(backgroundValues);
metrics.TargetStd = std(targetValues,0);
metrics.BackgroundStd = std(backgroundValues,0);
metrics.CNR = abs(metrics.TargetMean-metrics.BackgroundMean)/ ...
    sqrt(metrics.TargetStd^2+metrics.BackgroundStd^2+eps);
metrics.gCNR = generalizedCNR(targetValues,backgroundValues,numGCNRBins);
metrics.TargetSNR = metrics.TargetMean/(metrics.BackgroundStd+eps);
metrics.ContrastDB = 20*log10((metrics.TargetMean+eps)/ ...
    (metrics.BackgroundMean+eps));
metrics.ArtifactToTargetDB = 10*log10( ...
    (sum(backgroundValues.^2)+eps)/(sum(targetValues.^2)+eps));
end

%% ========================================================================
function value = generalizedCNR(targetValues, backgroundValues, numBins)
targetValues = targetValues(isfinite(targetValues));
backgroundValues = backgroundValues(isfinite(backgroundValues));
if isempty(targetValues) || isempty(backgroundValues)
    value = NaN;
    return;
end
minimum = min([targetValues(:);backgroundValues(:)]);
maximum = max([targetValues(:);backgroundValues(:)]);
if maximum-minimum <= eps
    value = 0;
    return;
end
edges = linspace(minimum,maximum,numBins+1);
pTarget = histcounts(targetValues,edges,'Normalization','probability');
pBackground = histcounts(backgroundValues,edges,'Normalization','probability');
value = 1-sum(min(pTarget,pBackground));
value = min(max(value,0),1);
end

%% ========================================================================
function normalized = normalizeAmplitude(image, roiMask)
image = abs(real(double(image)));
scale = max(image(roiMask));
if isempty(scale) || ~isfinite(scale) || scale <= eps
    normalized = zeros(size(image));
else
    normalized = min(max(image/scale,0),1);
end
normalized(~roiMask) = 0;
end

%% ========================================================================
function alpha = optimalGlobalScale(rawImage, gtP0, roiMask)
r = abs(real(double(rawImage)));
g = double(gtP0);
rv = r(roiMask);
gv = g(roiMask);
den = sum(rv.^2);
if isempty(den) || ~isfinite(den) || den <= eps
    alpha = NaN;
else
    alpha = sum(rv.*gv) / den;
end
end

%% ========================================================================
function metrics = computeAmplitudeFidelityMetrics(rawImage, gtP0, roiMask, ...
    targetMask, X, Y, bodyRadius, commonDASScale)
r = abs(real(double(rawImage)));
g = double(gtP0);

alpha = optimalGlobalScale(r, g, roiMask);
metrics.OptimalGlobalScale = alpha;
metrics.DASCommonScale = commonDASScale;
metrics.NRMSEOptimalScale = scaledNRMSE(r,g,roiMask,alpha);
metrics.TargetRecoveryOptimalScale = targetRecovery(r,g,targetMask,alpha);
metrics.NRMSEDASCommonScale = scaledNRMSE(r,g,roiMask,commonDASScale);
metrics.TargetRecoveryDASCommonScale = targetRecovery(r,g,targetMask,commonDASScale);

radius = hypot(X,Y);
edges = [0, bodyRadius/3, 2*bodyRadius/3, bodyRadius+eps];
metrics.TargetRecoveryInner = targetRecovery(r,g,targetMask & radius>=edges(1) & radius<edges(2),alpha);
metrics.TargetRecoveryMiddle = targetRecovery(r,g,targetMask & radius>=edges(2) & radius<edges(3),alpha);
metrics.TargetRecoveryOuter = targetRecovery(r,g,targetMask & radius>=edges(3) & radius<=edges(4),alpha);
end

%% ========================================================================
function value = scaledNRMSE(r,g,mask,alpha)
if ~isfinite(alpha) || ~any(mask(:))
    value = NaN;
    return;
end
rv = alpha*r(mask);
gv = g(mask);
value = norm(rv-gv) / max(norm(gv),eps);
end

%% ========================================================================
function value = targetRecovery(r,g,mask,alpha)
if ~isfinite(alpha) || nnz(mask)<1
    value = NaN;
    return;
end
den = mean(g(mask));
if ~isfinite(den) || abs(den)<=eps
    value = NaN;
else
    value = mean(alpha*r(mask))/den;
end
end

%% ========================================================================
function [hdMM,hd95MM] = boundaryDistancesMM(recMask,gtMask,pixelSize)
recEdge = bwperim(recMask,8);
gtEdge = bwperim(gtMask,8);
if ~any(recEdge(:)) && ~any(gtEdge(:))
    hdMM = 0;
    hd95MM = 0;
    return;
elseif ~any(recEdge(:)) || ~any(gtEdge(:))
    hdMM = Inf;
    hd95MM = Inf;
    return;
end
distanceToGT = bwdist(gtEdge);
distanceToRec = bwdist(recEdge);
d1 = distanceToGT(recEdge);
d2 = distanceToRec(gtEdge);
d = [d1(:);d2(:)];
hdMM = max(d)*pixelSize*1e3;
hd95MM = percentileLinear(d,95)*pixelSize*1e3;
end

%% ========================================================================
function value = percentileLinear(x,percentile)
x = sort(double(x(:)));
if isempty(x)
    value = NaN;
    return;
end
position = 1+(numel(x)-1)*percentile/100;
lo = floor(position);
hi = ceil(position);
if lo == hi
    value = x(lo);
else
    value = x(lo)+(position-lo)*(x(hi)-x(lo));
end
end

%% ========================================================================
function results = emptyMainResultsTable()
names = {'Phantom','SensorGeometry','Algorithm','RequestedSNR_dB', ...
    'MeasuredSNR_dB','NoiseRealization','KnownNoiseSigma_Pa', ...
    'EstimatedNoiseSigma_Pa','RMSE_ROI','PSNR_ROI_dB','SSIM','Dice', ...
    'Hausdorff_mm','HD95_mm','CNR','gCNR','TargetSNR','Contrast_dB', ...
    'ArtifactToTarget_dB','TargetMean','BackgroundMean','TargetStd', ...
    'BackgroundStd', ...
    'OptimalGlobalScale','NRMSE_OptimalScale','TargetRecovery_OptimalScale', ...
    'DASCommonScale','NRMSE_DASCommonScale','TargetRecovery_DASCommonScale', ...
    'TargetRecovery_Inner','TargetRecovery_Middle','TargetRecovery_Outer', ...
    'DelayRuntime_Sec','BeamformRuntime_Sec','TotalRuntime_Sec','NumSensors', ...
    'MVSubarrayLength','MVDiagonalLoading','SCFExponent','DASDSFExponent', ...
    'NCVDASExponent','NCVDASNoisePenalty','NoiseVarianceSource','PeakP0_Pa', ...
    'MeanFluence_Jm2','NumTargetPixels','NumBackgroundPixels','DelayModel','RunProfile'};
types = [repmat({'string'},1,3), repmat({'double'},1,39), ...
    {'string'}, repmat({'double'},1,4), {'string','string'}];
results = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function row = makeMainResultRow(phantomName, geometry, algorithm, ...
    requestedSNR, measuredSNR, realization, knownNoiseSigma, ...
    estimatedNoiseSigma, metrics, amplitudeMetrics, delayRuntime, beamformRuntime, ...
    totalRuntime, numSensors, cfg, gtP0, optical, targetMask, backgroundMask)
row = emptyMainResultsTable();
row(1,:) = {string(phantomName),string(geometry),string(algorithm), ...
    requestedSNR,measuredSNR,realization,knownNoiseSigma,estimatedNoiseSigma, ...
    metrics.RMSE,metrics.PSNR,metrics.SSIM,metrics.Dice,metrics.HausdorffMM, ...
    metrics.HD95MM,metrics.CNR,metrics.gCNR,metrics.TargetSNR, ...
    metrics.ContrastDB,metrics.ArtifactToTargetDB,metrics.TargetMean, ...
    metrics.BackgroundMean,metrics.TargetStd,metrics.BackgroundStd, ...
    amplitudeMetrics.OptimalGlobalScale,amplitudeMetrics.NRMSEOptimalScale, ...
    amplitudeMetrics.TargetRecoveryOptimalScale,amplitudeMetrics.DASCommonScale, ...
    amplitudeMetrics.NRMSEDASCommonScale,amplitudeMetrics.TargetRecoveryDASCommonScale, ...
    amplitudeMetrics.TargetRecoveryInner,amplitudeMetrics.TargetRecoveryMiddle, ...
    amplitudeMetrics.TargetRecoveryOuter, ...
    delayRuntime,beamformRuntime,totalRuntime,numSensors,cfg.mvSubarrayLength, ...
    cfg.mvDiagonalLoading,cfg.scfExponent,cfg.dasdsfExponent, ...
    cfg.ncVdasExponent,cfg.ncVdasNoisePenalty,"PretriggerMAD",max(gtP0(:)), ...
    mean(optical.fluenceJm2(optical.fluenceJm2>0)),nnz(targetMask), ...
    nnz(backgroundMask),string(cfg.delayModel),string(cfg.profile)};
end

%% ========================================================================
function tf = mainResultExists(results,phantom,geometry,snrDb,realization,algorithm)
if isempty(results)
    tf = false;
    return;
end
tf = any(results.Phantom == string(phantom) & ...
    results.SensorGeometry == string(geometry) & ...
    snrVectorEqual(results.RequestedSNR_dB,snrDb) & ...
    results.NoiseRealization == realization & ...
    results.Algorithm == string(algorithm));
end

%% ========================================================================
function tf = allMainResultsExist(results,phantom,geometry,snrDb,realization,algorithms)
tf = true;
for k = 1:numel(algorithms)
    if ~mainResultExists(results,phantom,geometry,snrDb,realization,algorithms(k))
        tf = false;
        return;
    end
end
end

%% ========================================================================
function mask = snrVectorEqual(values,target)
if isinf(target)
    mask = isinf(values) & sign(values)==sign(target);
else
    mask = abs(values-target)<1e-9;
end
end

%% ========================================================================
function tf = snrEqual(a,b)
if isinf(a) || isinf(b)
    tf = isinf(a) && isinf(b) && sign(a)==sign(b);
else
    tf = abs(a-b)<1e-9;
end
end

%% ========================================================================
function tf = anySNRMatch(value,array)
tf = false;
for k = 1:numel(array)
    if snrEqual(value,array(k))
        tf = true;
        return;
    end
end
end

%% ========================================================================
function text = formatSNR(value)
if isinf(value)
    text = 'Inf';
else
    text = sprintf('%g',value);
end
end

%% ========================================================================
function results = uniqueMainResults(results)
if isempty(results)
    return;
end
keys = results.Phantom+"|"+results.SensorGeometry+"|"+ ...
    string(results.RequestedSNR_dB)+"|"+string(results.NoiseRealization)+"|"+ ...
    results.Algorithm;
[~,uniqueIndex] = unique(keys,'stable');
results = results(uniqueIndex,:);
end

%% ========================================================================
function tableOut = emptyEquivalenceTable()
names = {'Phantom','SensorGeometry','RequestedSNR_dB','NoiseRealization', ...
    'RelativeL2Error','MaximumAbsoluteDifference','ImageCorrelation', ...
    'DirectRuntime_Sec','SMRuntime_Sec','Speedup_DirectOverSM', ...
    'MVSubarrayLength','MVDiagonalLoading'};
types = {'string','string','double','double','double','double','double', ...
    'double','double','double','double','double'};
tableOut = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function row = compareMVImages(phantom,geometry,snrDb,realization, ...
    directImage,smImage,directRuntime,smRuntime,L,loading)
relativeError = norm(directImage(:)-smImage(:),2)/ ...
    max(norm(directImage(:),2),eps);
maximumDifference = max(abs(directImage(:)-smImage(:)));
if std(directImage(:))<=eps || std(smImage(:))<=eps
    correlation = NaN;
else
    C = corrcoef(directImage(:),smImage(:));
    correlation = C(1,2);
end
row = emptyEquivalenceTable();
row(1,:) = {string(phantom),string(geometry),snrDb,realization, ...
    relativeError,maximumDifference,correlation,directRuntime,smRuntime, ...
    directRuntime/max(smRuntime,eps),L,loading};
end

%% ========================================================================
function tableOut = uniqueEquivalenceResults(tableOut)
if isempty(tableOut)
    return;
end
keys = tableOut.Phantom+"|"+tableOut.SensorGeometry+"|"+ ...
    string(tableOut.RequestedSNR_dB)+"|"+string(tableOut.NoiseRealization);
[~,idx] = unique(keys,'stable');
tableOut = tableOut(idx,:);
end

%% ========================================================================
function tableOut = emptySensitivityTable()
names = {'Phantom','SensorGeometry','RequestedSNR_dB','NoiseRealization', ...
    'ParameterFamily','Parameter1Name','Parameter1Value','Parameter2Name', ...
    'Parameter2Value','Algorithm','PSNR_ROI_dB','SSIM','Dice','HD95_mm', ...
    'CNR','gCNR','ArtifactToTarget_dB','Runtime_Sec'};
types = {'string','string','double','double','string','string','double', ...
    'string','double','string','double','double','double','double','double', ...
    'double','double','double'};
tableOut = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function results = runParameterSensitivity(delayed,gt,roiMask,targetMask, ...
    backgroundMask,context,cfg,phantomName,geometry,snrDb,realization)
fprintf('    Running non-iterative parameter-sensitivity experiment...\n');
results = emptySensitivityTable();

% SCF exponent
for value = cfg.sensitivitySCFExponent
    localCfg = cfg;
    localCfg.scfExponent = value;
    tic;
    vector = cheapBeamformerVector(double(delayed),"SCF-DAS", ...
        context.noiseVariance,localCfg);
    runtime = toc;
    results = [results; makeSensitivityRow(phantomName,geometry,snrDb, ...
        realization,"SCF_EXPONENT","p_SCF",value,"None",NaN, ...
        "SCF-DAS",vector,runtime,gt,roiMask,targetMask,backgroundMask,cfg)]; %#ok<AGROW>
end

% Bounded DASDSF exponent
for value = cfg.sensitivityDASDSFExponent
    localCfg = cfg;
    localCfg.dasdsfExponent = value;
    tic;
    vector = cheapBeamformerVector(double(delayed),"DASDSF", ...
        context.noiseVariance,localCfg);
    runtime = toc;
    results = [results; makeSensitivityRow(phantomName,geometry,snrDb, ...
        realization,"DASDSF_EXPONENT","q_DSF",value,"None",NaN, ...
        "DASDSF",vector,runtime,gt,roiMask,targetMask,backgroundMask,cfg)]; %#ok<AGROW>
end

% NC-VDAS exponent and noise penalty
for exponent = cfg.sensitivityNCVDASExponent
    for penalty = cfg.sensitivityNCVDASPenalty
        localCfg = cfg;
        localCfg.ncVdasExponent = exponent;
        localCfg.ncVdasNoisePenalty = penalty;
        tic;
        vector = cheapBeamformerVector(double(delayed),"NC-VDAS", ...
            context.noiseVariance,localCfg);
        runtime = toc;
        results = [results; makeSensitivityRow(phantomName,geometry,snrDb, ...
            realization,"NCVDAS_PARAMETERS","q_NCVDAS",exponent, ...
            "beta_noise",penalty,"NC-VDAS",vector,runtime,gt,roiMask, ...
            targetMask,backgroundMask,cfg)]; %#ok<AGROW>
    end
end

% MV subarray length and diagonal loading. For a circular full ring,
% cyclic subarray averaging makes the current SM-MV image algebraically
% collapse to DAS; therefore MV-parameter sensitivity is informative only
% for the non-circular limited-view geometry.
if ~context.isCircularArray
    for L = cfg.sensitivityMVSubarray
        if L > floor(size(delayed,1)/2)
            continue;
        end
        for loading = cfg.sensitivityMVLoading
            localCfg = cfg;
            localCfg.mvSubarrayLength = L;
            localCfg.mvDiagonalLoading = loading;
            tic;
            vector = spatiallySmoothedMVTiled(double(delayed),localCfg,context, ...
                "SHERMAN_MORRISON");
            runtime = toc;
            results = [results; makeSensitivityRow(phantomName,geometry,snrDb, ...
                realization,"MV_PARAMETERS","SubarrayLength",L, ...
                "DiagonalLoading",loading,"SM-MV",vector,runtime,gt,roiMask, ...
                targetMask,backgroundMask,cfg)]; %#ok<AGROW>
        end
    end
end
end

%% ========================================================================
function row = makeSensitivityRow(phantom,geometry,snrDb,realization, ...
    family,p1Name,p1Value,p2Name,p2Value,algorithm,vector,runtime,gt, ...
    roiMask,targetMask,backgroundMask,cfg)
rec = normalizeAmplitude(reshape(vector,cfg.Nx,cfg.Ny),roiMask);
metrics = computeArticleMetrics(rec,gt,roiMask,targetMask,backgroundMask, ...
    cfg.segmentationThreshold,cfg.dx,cfg.gCNRBins);
row = emptySensitivityTable();
row(1,:) = {string(phantom),string(geometry),snrDb,realization, ...
    string(family),string(p1Name),p1Value,string(p2Name),p2Value, ...
    string(algorithm),metrics.PSNR,metrics.SSIM,metrics.Dice,metrics.HD95MM, ...
    metrics.CNR,metrics.gCNR,metrics.ArtifactToTargetDB,runtime};
end

%% ========================================================================
function tableOut = uniqueSensitivityResults(tableOut)
if isempty(tableOut)
    return;
end
keys = tableOut.Phantom+"|"+tableOut.SensorGeometry+"|"+ ...
    string(tableOut.RequestedSNR_dB)+"|"+string(tableOut.NoiseRealization)+"|"+ ...
    tableOut.ParameterFamily+"|"+string(tableOut.Parameter1Value)+"|"+ ...
    string(tableOut.Parameter2Value)+"|"+tableOut.Algorithm;
[~,idx] = unique(keys,'stable');
tableOut = tableOut(idx,:);
end

%% ========================================================================
function summary = makeArticleSummary(results,cfg)
if isempty(results)
    summary = table();
    return;
end
cohortNames = ["CHEAP_ALL_PHANTOMS","ADAPTIVE_COMMON"];
records = cell(0,1);
row = 1;
metricNames = ["PSNR_ROI_dB","SSIM","Dice","HD95_mm","CNR","gCNR", ...
    "ArtifactToTarget_dB","NRMSE_OptimalScale","TargetRecovery_OptimalScale", ...
    "NRMSE_DASCommonScale","TargetRecovery_DASCommonScale", ...
    "TargetRecovery_Inner","TargetRecovery_Middle","TargetRecovery_Outer", ...
    "TotalRuntime_Sec"];
for cohortIndex = 1:numel(cohortNames)
    cohort = cohortNames(cohortIndex);
    switch cohort
        case "CHEAP_ALL_PHANTOMS"
            cohortRows = ismember(results.Algorithm,cfg.cheapAlgorithms);
        case "ADAPTIVE_COMMON"
            scheduleRows = false(height(results),1);
            for rowIndex = 1:height(results)
                scheduleRows(rowIndex) = ...
                    anySNRMatch(results.RequestedSNR_dB(rowIndex),cfg.adaptiveSNRdB) && ...
                    (isinf(results.RequestedSNR_dB(rowIndex)) || ...
                    results.NoiseRealization(rowIndex)<=cfg.adaptiveNoiseRealizations);
            end
            cohortRows = ismember(results.Phantom,cfg.adaptivePhantomNames) & ...
                results.Algorithm~="SS-MV" & scheduleRows;
            degenerateFullRingMV = results.SensorGeometry=="CIRCULAR_FULL" & ...
                contains(results.Algorithm,"SM-MV");
            cohortRows = cohortRows & ~degenerateFullRingMV;
        otherwise
            cohortRows = true(height(results),1);
    end
    cohortResults = results(cohortRows,:);
    if isempty(cohortResults)
        continue;
    end
    algorithms = unique(cohortResults.Algorithm,'stable');
    geometries = unique(cohortResults.SensorGeometry,'stable');
    snrValues = unique(cohortResults.RequestedSNR_dB,'stable');
    for g = 1:numel(geometries)
        for s = 1:numel(snrValues)
            for a = 1:numel(algorithms)
                idx = cohortResults.SensorGeometry==geometries(g) & ...
                    snrVectorEqual(cohortResults.RequestedSNR_dB,snrValues(s)) & ...
                    cohortResults.Algorithm==algorithms(a);
                if ~any(idx)
                    continue;
                end
                phantomLabels = cohortResults.Phantom(idx);
                uniquePhantoms = unique(phantomLabels,'stable');
                rec = struct();
                rec.EvaluationCohort = cohort;
                rec.SummaryUnit = "PHANTOM_CLUSTER";
                rec.SensorGeometry = geometries(g);
                rec.RequestedSNR_dB = snrValues(s);
                rec.Algorithm = algorithms(a);
                rec.PhantomCount = numel(uniquePhantoms);
                rec.SampleCount = nnz(idx);
                for metricIndex = 1:numel(metricNames)
                    values = cohortResults.(metricNames(metricIndex))(idx);
                    clusterValues = zeros(numel(uniquePhantoms),1);
                    for clusterIndex = 1:numel(uniquePhantoms)
                        clusterRows = phantomLabels==uniquePhantoms(clusterIndex);
                        clusterValues(clusterIndex) = mean(values(clusterRows),'omitnan');
                    end
                    prefix = metricNames(metricIndex);
                    rec.(prefix+"_Mean") = mean(clusterValues,'omitnan');
                    rec.(prefix+"_Std") = std(clusterValues,0,'omitnan');
                    rec.(prefix+"_Median") = median(clusterValues,'omitnan');
                    rec.(prefix+"_IQR") = percentileLinear(clusterValues,75)- ...
                        percentileLinear(clusterValues,25);
                    n = nnz(isfinite(clusterValues));
                    halfWidth = 1.96*std(clusterValues,0,'omitnan')/sqrt(max(n,1));
                    rec.(prefix+"_CI95Low") = rec.(prefix+"_Mean")-halfWidth;
                    rec.(prefix+"_CI95High") = rec.(prefix+"_Mean")+halfWidth;
                end
                records{row,1} = rec; %#ok<AGROW>
                row = row+1;
            end
        end
    end
end
if isempty(records)
    summary = table();
else
    summary = struct2table(vertcat(records{:}));
end
end

%% ========================================================================
function paired = makePairedComparisons(results,cfg)
metricNames = ["PSNR_ROI_dB","SSIM","Dice","HD95_mm","CNR","gCNR", ...
    "ArtifactToTarget_dB","NRMSE_OptimalScale","NRMSE_DASCommonScale"];
higherIsBetter = [true,true,true,false,true,true,false,false,false];
algorithms = setdiff(unique(results.Algorithm,'stable'),["DAS","SS-MV"],'stable');
geometries = unique(results.SensorGeometry,'stable');
snrValues = unique(results.RequestedSNR_dB,'stable');
records = cell(0,1);
row = 1;
rng(cfg.statisticsSeed,'twister');
for g = 1:numel(geometries)
    for s = 1:numel(snrValues)
        baselineRows = results.SensorGeometry==geometries(g) & ...
            snrVectorEqual(results.RequestedSNR_dB,snrValues(s)) & ...
            results.Algorithm=="DAS";
        baseline = results(baselineRows,:);
        if isempty(baseline)
            continue;
        end
        baselineKey = makePairKey(baseline);
        for a = 1:numel(algorithms)
            candidateRows = results.SensorGeometry==geometries(g) & ...
                snrVectorEqual(results.RequestedSNR_dB,snrValues(s)) & ...
                results.Algorithm==algorithms(a);
            candidate = results(candidateRows,:);
            if isempty(candidate)
                continue;
            end
            candidateKey = makePairKey(candidate);
            [~,ib,ic] = intersect(baselineKey,candidateKey,'stable');
            if numel(ib)<2
                continue;
            end
            pairedPhantoms = candidate.Phantom(ic);
            for metricIndex = 1:numel(metricNames)
                baselineValues = baseline.(metricNames(metricIndex))(ib);
                candidateValues = candidate.(metricNames(metricIndex))(ic);
                valid = isfinite(baselineValues) & isfinite(candidateValues);
                baselineValues = baselineValues(valid);
                candidateValues = candidateValues(valid);
                phantomLabels = pairedPhantoms(valid);
                if numel(baselineValues)<2
                    continue;
                end

                % Average repeated noise realizations within each phantom.
                uniquePhantoms = unique(phantomLabels,'stable');
                clusterBaseline = zeros(numel(uniquePhantoms),1);
                clusterCandidate = zeros(numel(uniquePhantoms),1);
                for clusterIndex = 1:numel(uniquePhantoms)
                    clusterRows = phantomLabels==uniquePhantoms(clusterIndex);
                    clusterBaseline(clusterIndex) = mean(baselineValues(clusterRows));
                    clusterCandidate(clusterIndex) = mean(candidateValues(clusterRows));
                end
                clusterRawDifference = clusterCandidate-clusterBaseline;
                if higherIsBetter(metricIndex)
                    clusterOrientedDifference = clusterRawDifference;
                else
                    clusterOrientedDifference = -clusterRawDifference;
                end

                [ciLow,ciHigh] = bootstrapMeanCI(clusterRawDifference, ...
                    cfg.bootstrapCount);
                pValue = pairedSignFlipP(clusterOrientedDifference, ...
                    cfg.permutationCount);
                effectDz = mean(clusterOrientedDifference)/ ...
                    max(std(clusterOrientedDifference,0),eps);

                rec = struct();
                rec.EvaluationUnit = "PHANTOM_CLUSTER";
                rec.SensorGeometry = geometries(g);
                rec.RequestedSNR_dB = snrValues(s);
                rec.Algorithm = algorithms(a);
                rec.Baseline = "DAS";
                rec.Metric = metricNames(metricIndex);
                rec.NumMatchedRealizations = numel(baselineValues);
                rec.NumPhantoms = numel(uniquePhantoms);
                rec.MeanCandidate = mean(clusterCandidate);
                rec.MeanBaseline = mean(clusterBaseline);
                rec.MeanDifference_CandidateMinusBaseline = ...
                    mean(clusterRawDifference);
                rec.CI95Low = ciLow;
                rec.CI95High = ciHigh;
                rec.OrientedEffectSizeDz = effectDz;
                rec.PermutationP = pValue;
                rec.HigherIsBetter = higherIsBetter(metricIndex);
                records{row,1} = rec; %#ok<AGROW>
                row = row+1;
            end
        end
    end
end
if isempty(records)
    paired = table();
else
    paired = struct2table(vertcat(records{:}));
    paired.FDR_Q = benjaminiHochberg(paired.PermutationP);
end
end

%% ========================================================================
function key = makePairKey(tableIn)
key = tableIn.Phantom+"|"+tableIn.SensorGeometry+"|"+ ...
    string(tableIn.RequestedSNR_dB)+"|"+string(tableIn.NoiseRealization);
end

%% ========================================================================
function [low,high] = bootstrapMeanCI(differences,count)
differences = differences(:);
n = numel(differences);
if n<2 || count<10
    low = NaN;
    high = NaN;
    return;
end
means = zeros(count,1);
chunk = 500;
filled = 0;
while filled<count
    current = min(chunk,count-filled);
    indices = randi(n,n,current);
    means(filled+(1:current)) = mean(differences(indices),1).';
    filled = filled+current;
end
low = percentileLinear(means,2.5);
high = percentileLinear(means,97.5);
end

%% ========================================================================
function pValue = pairedSignFlipP(orientedDifferences,count)
d = orientedDifferences(:);
d = d(isfinite(d));
if isempty(d)
    pValue = NaN;
    return;
end
observed = abs(mean(d));
n = numel(d);
if n<=16
    % Exact paired randomization test at the phantom-cluster level.
    numberPatterns = 2^n;
    integers = uint32((0:numberPatterns-1)');
    signs = ones(n,numberPatterns);
    for bitIndex = 1:n
        bitValues = bitget(integers,bitIndex).';
        signs(bitIndex,:) = 2*double(bitValues)-1;
    end
    permutedMeans = abs(mean(d.*signs,1));
    pValue = nnz(permutedMeans>=observed-10*eps)/numberPatterns;
else
    extreme = 0;
    chunk = 500;
    completed = 0;
    while completed<count
        current = min(chunk,count-completed);
        signs = 2*(rand(n,current)>0.5)-1;
        permutedMeans = abs(mean(d.*signs,1));
        extreme = extreme+nnz(permutedMeans>=observed-10*eps);
        completed = completed+current;
    end
    pValue = (extreme+1)/(count+1);
end
end

%% ========================================================================
function qValues = benjaminiHochberg(pValues)
p = pValues(:);
qValues = NaN(size(p));
valid = find(isfinite(p));
if isempty(valid)
    return;
end
[pSorted,order] = sort(p(valid));
m = numel(pSorted);
qSorted = pSorted.*m./(1:m)';
for k = m-1:-1:1
    qSorted(k) = min(qSorted(k),qSorted(k+1));
end
qSorted = min(qSorted,1);
temp = zeros(m,1);
temp(order) = qSorted;
qValues(valid) = temp;
end

%% ========================================================================
function rankings = makeAlgorithmRankings(summary)
if isempty(summary)
    rankings = table();
    return;
end
if any(summary.EvaluationCohort=="ADAPTIVE_COMMON")
    summary = summary(summary.EvaluationCohort=="ADAPTIVE_COMMON",:);
else
    summary = summary(summary.EvaluationCohort==summary.EvaluationCohort(1),:);
end
geometries = unique(summary.SensorGeometry,'stable');
snrValues = unique(summary.RequestedSNR_dB,'stable');
records = cell(0,1);
row = 1;
for g = 1:numel(geometries)
    for s = 1:numel(snrValues)
        idx = summary.SensorGeometry==geometries(g) & ...
            snrVectorEqual(summary.RequestedSNR_dB,snrValues(s));
        block = summary(idx,:);
        if isempty(block)
            continue;
        end
        maximumPhantomCount = max(block.PhantomCount);
        block = block(block.PhantomCount==maximumPhantomCount,:);
        rankPSNR = tiedRanks(-block.PSNR_ROI_dB_Mean);
        rankSSIM = tiedRanks(-block.SSIM_Mean);
        rankGCNR = tiedRanks(-block.gCNR_Mean);
        rankHD95 = tiedRanks(block.HD95_mm_Mean);
        aggregate = mean([rankPSNR,rankSSIM,rankGCNR,rankHD95],2);
        for k = 1:height(block)
            rec = struct();
            rec.EvaluationCohort = block.EvaluationCohort(k);
            rec.SensorGeometry = geometries(g);
            rec.RequestedSNR_dB = snrValues(s);
            rec.Algorithm = block.Algorithm(k);
            rec.RankPSNR = rankPSNR(k);
            rec.RankSSIM = rankSSIM(k);
            rec.RankGCNR = rankGCNR(k);
            rec.RankHD95 = rankHD95(k);
            rec.AggregateRank = aggregate(k);
            records{row,1} = rec; %#ok<AGROW>
            row = row+1;
        end
    end
end
if isempty(records)
    rankings = table();
else
    rankings = struct2table(vertcat(records{:}));
    rankings = sortrows(rankings,{'SensorGeometry','RequestedSNR_dB','AggregateRank'});
end
end

%% ========================================================================
function r = tiedRanks(values)
values = values(:);
[sorted,order] = sort(values);
rSorted = zeros(size(sorted));
i = 1;
while i<=numel(sorted)
    j = i;
    while j<numel(sorted) && abs(sorted(j+1)-sorted(i))<1e-12
        j = j+1;
    end
    rSorted(i:j) = mean(i:j);
    i = j+1;
end
r = zeros(size(values));
r(order) = rSorted;
end

%% ========================================================================
function best = makeBestAlgorithmTable(summary)
if isempty(summary)
    best = table();
    return;
end
if any(summary.EvaluationCohort=="ADAPTIVE_COMMON")
    summary = summary(summary.EvaluationCohort=="ADAPTIVE_COMMON",:);
else
    summary = summary(summary.EvaluationCohort==summary.EvaluationCohort(1),:);
end
geometries = unique(summary.SensorGeometry,'stable');
snrValues = unique(summary.RequestedSNR_dB,'stable');
metrics = ["PSNR_ROI_dB_Mean","SSIM_Mean","gCNR_Mean","HD95_mm_Mean", ...
    "TotalRuntime_Sec_Median"];
higher = [true,true,true,false,false];
records = cell(0,1);
row = 1;
for g = 1:numel(geometries)
    for s = 1:numel(snrValues)
        block = summary(summary.SensorGeometry==geometries(g) & ...
            snrVectorEqual(summary.RequestedSNR_dB,snrValues(s)),:);
        maximumPhantomCount = max(block.PhantomCount);
        block = block(block.PhantomCount==maximumPhantomCount,:);
        for m = 1:numel(metrics)
            values = block.(metrics(m));
            if higher(m)
                [bestValue,index] = max(values);
            else
                [bestValue,index] = min(values);
            end
            rec = struct();
            rec.EvaluationCohort = block.EvaluationCohort(index);
            rec.SensorGeometry = geometries(g);
            rec.RequestedSNR_dB = snrValues(s);
            rec.Metric = metrics(m);
            rec.BestAlgorithm = block.Algorithm(index);
            rec.BestValue = bestValue;
            records{row,1} = rec; %#ok<AGROW>
            row = row+1;
        end
    end
end
if isempty(records)
    best = table();
else
    best = struct2table(vertcat(records{:}));
end
end

%% ========================================================================
function results = runResolutionBenchmark(background,kgrid,karray, ...
    travelTimesFull,limitedIndices,X,Y,xVec,yVec,tArray,roiMask,cfg)
outputDir = fullfile(cfg.outputDir,'resolution_benchmark');
ensureDirectory(outputDir);
checkpointFile = fullfile(outputDir,'resolution_checkpoint.mat');
if cfg.resume && isfile(checkpointFile)
    loadedCheckpoint = load(checkpointFile,'results');
    if isfield(loadedCheckpoint,'results')
        results = loadedCheckpoint.results;
    else
        results = emptyResolutionTable();
    end
else
    results = emptyResolutionTable();
end

for pointIndex = 1:size(cfg.resolutionPointLocationsMM,1)
    pointMM = cfg.resolutionPointLocationsMM(pointIndex,:);
    pointM = pointMM*1e-3;
    sigma = cfg.resolutionPointSigmaMM*1e-3;
    p0 = exp(-((X-pointM(1)).^2+(Y-pointM(2)).^2)/(2*sigma^2));
    p0 = 1000*p0/max(p0(:));

    dataFile = fullfile(outputDir,sprintf('point_%02d_clean_data.mat',pointIndex));
    if cfg.resume && isfile(dataFile)
        loaded = load(dataFile,'cleanFullData');
        cleanFullData = loaded.cleanFullData;
    else
        cleanFullData = single(simulateKWaveFiniteArray(p0,background, ...
            kgrid,karray,cfg));
        save(dataFile,'cleanFullData','p0','pointMM','-v7.3');
    end

    for geometry = cfg.sensorGeometries
        [sensorIndices,isCircular] = geometryIndices(geometry,cfg,limitedIndices);
        cleanData = double(cleanFullData(sensorIndices,:));
        travelTimes = travelTimesFull(sensorIndices,:);

        for snrDb = cfg.resolutionSNRdB
            if isinf(snrDb)
                repetitions = 1;
            else
                repetitions = cfg.resolutionRealizations;
            end
            for realization = 1:repetitions
                rng(caseSeed(cfg.masterSeed+700001,pointIndex,geometry,snrDb,realization), ...
                    'twister');
                [noisyData,~,measuredSNR] = addNoiseAtRequestedSNR(cleanData,snrDb);
                noiseVar = estimatePretriggerNoiseVariance(noisyData, ...
                    cfg.noiseEstimationSamples);
                tic;
                delayed = sampleDelayedData(noisyData,tArray,travelTimes);
                delayRuntime = toc;
                context.isCircularArray = isCircular;
                context.pixelMask = roiMask;
                context.noiseVariance = noiseVar;
                context.numSensors = numel(sensorIndices);
                outputs = reconstructAlgorithmSet(delayed,cfg.resolutionAlgorithms, ...
                    cfg,context);

                for a = 1:numel(outputs)
                    algorithmName = string(outputs(a).Algorithm);
                    if resolutionResultExistsV91(results, pointIndex, geometry, ...
                            algorithmName, snrDb, realization)
                        continue;
                    end
                    rec = normalizeAmplitude(outputs(a).RawImage,roiMask);
                    psf = measurePointSpread(rec,xVec,yVec,pointM,cfg);
                    row = emptyResolutionTable();
                    row(1,:) = {pointIndex,pointMM(1),pointMM(2),string(geometry), ...
                        string(outputs(a).Algorithm),snrDb,measuredSNR,realization, ...
                        psf.FWHMRadialMM,psf.FWHMTangentialMM, ...
                        psf.PSLRadialDB,psf.PSLTangentialDB, ...
                        psf.PeakLocalizationErrorMM,delayRuntime, ...
                        outputs(a).BeamformRuntimeSec, ...
                        delayRuntime+outputs(a).BeamformRuntimeSec};
                    results = [results;row]; %#ok<AGROW>

                    if realization==1
                        savePSFProfileFigure(psf,rec,xVec,yVec,pointIndex,pointMM, ...
                            geometry,string(outputs(a).Algorithm),snrDb,outputDir);
                    end
                end

                % Preserve partial resolution results after each realization.
                save(fullfile(outputDir,'resolution_checkpoint.mat'), ...
                    'results','-v7.3');
                writetable(results, ...
                    fullfile(outputDir,'resolution_results_checkpoint.csv'));
            end
        end
    end
end
end

%% ========================================================================
function tf = resolutionResultExistsV91(results, pointIndex, geometry, ...
    algorithm, snrDb, realization)
if isempty(results)
    tf = false;
    return;
end
if isinf(snrDb)
    snrMask = isinf(results.RequestedSNR_dB) & ...
        sign(results.RequestedSNR_dB) == sign(snrDb);
else
    snrMask = abs(results.RequestedSNR_dB-snrDb) < 1e-9;
end
tf = any(results.PointIndex == pointIndex & ...
    results.SensorGeometry == string(geometry) & ...
    results.Algorithm == string(algorithm) & snrMask & ...
    results.NoiseRealization == realization);
end

%% ========================================================================
function tableOut = emptyResolutionTable()
names = {'PointIndex','PointX_mm','PointY_mm','SensorGeometry','Algorithm', ...
    'RequestedSNR_dB','MeasuredSNR_dB','NoiseRealization', ...
    'FWHM_Radial_mm','FWHM_Tangential_mm','PSL_Radial_dB', ...
    'PSL_Tangential_dB','PeakLocalizationError_mm','DelayRuntime_Sec', ...
    'BeamformRuntime_Sec','TotalRuntime_Sec'};
types = {'double','double','double','string','string','double','double', ...
    'double','double','double','double','double','double','double','double','double'};
tableOut = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function psf = measurePointSpread(rec,xVec,yVec,pointM,cfg)
if hypot(pointM(1),pointM(2))<1e-12
    radialUnit = [1,0];
else
    radialUnit = pointM/hypot(pointM(1),pointM(2));
end
tangentialUnit = [-radialUnit(2),radialUnit(1)];
s = linspace(-cfg.profileHalfWidthMM,cfg.profileHalfWidthMM, ...
    cfg.profileSamples)*1e-3;

radialX = pointM(1)+s*radialUnit(1);
radialY = pointM(2)+s*radialUnit(2);
tangentialX = pointM(1)+s*tangentialUnit(1);
tangentialY = pointM(2)+s*tangentialUnit(2);

radialProfile = interp2(yVec,xVec,rec,radialY,radialX,'linear',0);
tangentialProfile = interp2(yVec,xVec,rec,tangentialY,tangentialX,'linear',0);
radialProfile = radialProfile/max(radialProfile+eps);
tangentialProfile = tangentialProfile/max(tangentialProfile+eps);

psf.FWHMRadialMM = profileFWHM(s*1e3,radialProfile);
psf.FWHMTangentialMM = profileFWHM(s*1e3,tangentialProfile);
psf.PSLRadialDB = profilePSL(s*1e3,radialProfile,psf.FWHMRadialMM);
psf.PSLTangentialDB = profilePSL(s*1e3,tangentialProfile,psf.FWHMTangentialMM);

[~,linearIndex] = max(rec(:));
[row,col] = ind2sub(size(rec),linearIndex);
peakX = xVec(row);
peakY = yVec(col);
psf.PeakLocalizationErrorMM = hypot(peakX-pointM(1),peakY-pointM(2))*1e3;
psf.DistanceMM = s*1e3;
psf.RadialProfile = radialProfile;
psf.TangentialProfile = tangentialProfile;
end

%% ========================================================================
function width = profileFWHM(distanceMM,profile)
profile = double(profile(:));
distanceMM = double(distanceMM(:));
[peak,index] = max(profile);
if peak<=0
    width = NaN;
    return;
end
half = 0.5*peak;
leftIndex = find(profile(1:index)<=half,1,'last');
rightRelative = find(profile(index:end)<=half,1,'first');
if isempty(leftIndex) || isempty(rightRelative)
    width = NaN;
    return;
end
rightIndex = index+rightRelative-1;
leftCross = interpolateCrossing(distanceMM(leftIndex:leftIndex+1), ...
    profile(leftIndex:leftIndex+1),half);
rightCross = interpolateCrossing(distanceMM(rightIndex-1:rightIndex), ...
    profile(rightIndex-1:rightIndex),half);
width = rightCross-leftCross;
end

%% ========================================================================
function xCross = interpolateCrossing(x,y,target)
if abs(y(2)-y(1))<=eps
    xCross = mean(x);
else
    xCross = x(1)+(target-y(1))*(x(2)-x(1))/(y(2)-y(1));
end
end

%% ========================================================================
function pslDb = profilePSL(distanceMM,profile,fwhmMM)
profile = abs(double(profile(:)));
distanceMM = double(distanceMM(:));
[peak,index] = max(profile); %#ok<ASGLU>
if peak<=0 || ~isfinite(fwhmMM)
    pslDb = NaN;
    return;
end
peakDistance = distanceMM(index);
exclude = abs(distanceMM-peakDistance)<=0.75*fwhmMM;
sideValues = profile(~exclude);
if isempty(sideValues)
    pslDb = NaN;
else
    pslDb = 20*log10(max(sideValues)/peak+eps);
end
end

%% ========================================================================
function results = runDelayMismatchBenchmark(background, phantomContrasts, ...
    phantomNames, selectedPhantomIndices, sensorXFull, sensorYFull, ...
    limitedIndices, X, Y, xVec, yVec, tArray, roiMask, cfg)
outputDir = fullfile(cfg.outputDir,'reviewer_delay_mismatch');
ensureDirectory(outputDir);
results = emptyDelayMismatchTable();

selectedNames = string(phantomNames(selectedPhantomIndices));
for requestedName = cfg.mismatchPhantomNames
    localPos = find(selectedNames==requestedName,1);
    if isempty(localPos)
        warning('PACT:MismatchPhantomUnavailable', ...
            'Mismatch phantom %s was not part of the selected main run.',char(requestedName));
        continue;
    end
    p = selectedPhantomIndices(localPos);
    phantomName = phantomNames(p);
    contrast = phantomContrasts{p};
    optical = buildOpticalAndPressureMaps(contrast,background,X,Y,cfg);
    gtP0 = optical.p0Pa;
    gtNorm = normalizeAmplitude(gtP0,roiMask);
    [targetMask,backgroundMask] = makeEvaluationMasks(gtNorm,roiMask, ...
        cfg.segmentationThreshold,cfg.backgroundGtMaximum, ...
        cfg.backgroundMarginMM,cfg.dx);

    dataFile = fullfile(cfg.outputDir,char(phantomName), ...
        'cached_forward_data','clean_full_ring_data.mat');
    if ~isfile(dataFile)
        warning('PACT:MissingCachedData','Missing cached data: %s',dataFile);
        continue;
    end
    loaded = load(dataFile,'cleanFullData');
    cleanFullData = double(loaded.cleanFullData);

    for condition = cfg.mismatchConditions
        [assumedMap,delayMode,homogeneousC] = makeMismatchSoundSpeed( ...
            background.soundSpeed,X,Y,xVec,yVec,condition,cfg);
        localCfg = cfg;
        localCfg.delayModel = delayMode;
        localCfg.homogeneousReconSoundSpeed = homogeneousC;
        travelTimesAll = computeTravelTimes(X,Y,xVec,yVec,sensorXFull, ...
            sensorYFull,assumedMap,localCfg);

        for geometry = cfg.mismatchGeometries
            [sensorIndices,isCircular] = geometryIndices(geometry,cfg,limitedIndices);
            cleanData = cleanFullData(sensorIndices,:);
            travelTimes = travelTimesAll(sensorIndices,:);
            for snrDb = cfg.mismatchSNRdB
                reps = cfg.mismatchRealizations;
                if isinf(snrDb), reps = 1; end
                for realization = 1:reps
                    rng(caseSeed(cfg.masterSeed+810001,p,geometry,snrDb,realization), ...
                        'twister');
                    [noisyData,~,measuredSNR] = addNoiseAtRequestedSNR(cleanData,snrDb);
                    noiseVar = estimatePretriggerNoiseVariance(noisyData, ...
                        cfg.noiseEstimationSamples);
                    delayed = sampleDelayedData(noisyData,tArray,travelTimes);
                    context.isCircularArray = isCircular;
                    context.pixelMask = roiMask;
                    context.noiseVariance = noiseVar;
                    context.numSensors = numel(sensorIndices);
                    algorithms = cfg.mismatchAlgorithms;
                    % Full-ring SM-MV is structurally degenerate to DAS for the
                    % current cyclic subarray averaging; retain only DAS/direct
                    % nonlinear methods in the full-ring mismatch ranking.
                    if isCircular
                        algorithms = algorithms(~contains(algorithms,"SM-MV"));
                    end
                    outputs = reconstructAlgorithmSet(delayed,algorithms,cfg,context);
                    for a = 1:numel(outputs)
                        recNorm = normalizeAmplitude(outputs(a).RawImage,roiMask);
                        metrics = computeArticleMetrics(recNorm,gtNorm,roiMask, ...
                            targetMask,backgroundMask,cfg.segmentationThreshold, ...
                            cfg.dx,cfg.gCNRBins);
                        row = emptyDelayMismatchTable();
                        row(1,:) = {string(phantomName),string(condition), ...
                            string(geometry),string(outputs(a).Algorithm),snrDb, ...
                            measuredSNR,realization,metrics.PSNR,metrics.SSIM, ...
                            metrics.Dice,metrics.HD95MM,metrics.gCNR, ...
                            metrics.ArtifactToTargetDB};
                        results = [results;row]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
end

%% ========================================================================
function tableOut = emptyDelayMismatchTable()
names = {'Phantom','MismatchCondition','SensorGeometry','Algorithm', ...
    'RequestedSNR_dB','MeasuredSNR_dB','NoiseRealization','PSNR_ROI_dB', ...
    'SSIM','Dice','HD95_mm','gCNR','ArtifactToTarget_dB'};
types = [{'string','string','string','string'},repmat({'double'},1,9)];
tableOut = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function [assumedMap,delayMode,homogeneousC] = makeMismatchSoundSpeed( ...
    trueMap,X,Y,xVec,yVec,condition,cfg)
condition = upper(string(condition));
homogeneousC = cfg.homogeneousReconSoundSpeed;
assumedMap = double(trueMap);
delayMode = "STRAIGHT_RAY_HETEROGENEOUS";

switch condition
    case "BASELINE_HETERO"
        % unchanged
    case "HOMOGENEOUS_1540"
        delayMode = "HOMOGENEOUS";
        homogeneousC = 1540;
        assumedMap(:) = homogeneousC;
    case "GLOBAL_MINUS_2PCT"
        assumedMap = (1-cfg.mismatchGlobalFraction)*double(trueMap);
    case "GLOBAL_PLUS_2PCT"
        assumedMap = (1+cfg.mismatchGlobalFraction)*double(trueMap);
    case "MISREG_X_0P5MM"
        shift = cfg.mismatchShiftMM*1e-3;
        assumedMap = interp2(yVec,xVec,double(trueMap),Y,X-shift, ...
            'linear',cfg.homogeneousReconSoundSpeed);
    case "LOWRES_X4"
        factor = max(2,round(cfg.mismatchDownsampleFactor));
        small = imresize(double(trueMap),1/factor,'bilinear');
        assumedMap = imresize(small,size(trueMap),'bilinear');
    case "UNDERSEGMENT_2PX"
        c0 = cfg.homogeneousReconSoundSpeed;
        heteroMask = abs(double(trueMap)-c0) > 5;
        radius = max(1,round(cfg.mismatchSegmentationPixels));
        eroded = imerode(heteroMask,strel('disk',radius,0));
        assumedMap = c0*ones(size(trueMap));
        assumedMap(eroded) = double(trueMap(eroded));
    otherwise
        error('Unknown mismatch condition: %s',char(condition));
end
assumedMap = max(1300,min(1800,assumedMap));
end

%% ========================================================================
function results = runNoiseRobustnessBenchmark(background, phantomContrasts, ...
    phantomNames, selectedPhantomIndices, limitedIndices, X, Y, tArray, ...
    roiMask, travelTimesFull, cfg)
outputDir = fullfile(cfg.outputDir,'reviewer_noise_robustness');
ensureDirectory(outputDir);
results = emptyNoiseRobustnessTable();
selectedNames = string(phantomNames(selectedPhantomIndices));

for requestedName = cfg.noiseRobustnessPhantomNames
    localPos = find(selectedNames==requestedName,1);
    if isempty(localPos), continue; end
    p = selectedPhantomIndices(localPos);
    phantomName = phantomNames(p);
    contrast = phantomContrasts{p};
    optical = buildOpticalAndPressureMaps(contrast,background,X,Y,cfg);
    gtNorm = normalizeAmplitude(optical.p0Pa,roiMask);
    [targetMask,backgroundMask] = makeEvaluationMasks(gtNorm,roiMask, ...
        cfg.segmentationThreshold,cfg.backgroundGtMaximum, ...
        cfg.backgroundMarginMM,cfg.dx);

    dataFile = fullfile(cfg.outputDir,char(phantomName), ...
        'cached_forward_data','clean_full_ring_data.mat');
    if ~isfile(dataFile), continue; end
    loaded = load(dataFile,'cleanFullData');
    cleanFullData = double(loaded.cleanFullData);

    for geometry = cfg.noiseRobustnessGeometries
        [sensorIndices,isCircular] = geometryIndices(geometry,cfg,limitedIndices);
        cleanData = cleanFullData(sensorIndices,:);
        travelTimes = travelTimesFull(sensorIndices,:);
        for noiseModel = cfg.noiseRobustnessModels
            for realization = 1:cfg.noiseRobustnessRealizations
                rng(caseSeed(cfg.masterSeed+820001,p,geometry, ...
                    cfg.noiseRobustnessSNRdB,realization)+sum(double(char(noiseModel))), ...
                    'twister');
                [noisyData,noiseSigma,measuredSNR] = addReviewerNoise( ...
                    cleanData,cfg.noiseRobustnessSNRdB,noiseModel,cfg);
                noiseVar = estimatePretriggerNoiseVariance(noisyData, ...
                    cfg.noiseEstimationSamples);
                delayed = sampleDelayedData(noisyData,tArray,travelTimes);
                context.isCircularArray = isCircular;
                context.pixelMask = roiMask;
                context.noiseVariance = noiseVar;
                context.numSensors = numel(sensorIndices);
                algorithms = cfg.noiseRobustnessAlgorithms;
                if isCircular
                    algorithms = algorithms(~contains(algorithms,"SM-MV"));
                end
                outputs = reconstructAlgorithmSet(delayed,algorithms,cfg,context);
                for a = 1:numel(outputs)
                    recNorm = normalizeAmplitude(outputs(a).RawImage,roiMask);
                    metrics = computeArticleMetrics(recNorm,gtNorm,roiMask, ...
                        targetMask,backgroundMask,cfg.segmentationThreshold, ...
                        cfg.dx,cfg.gCNRBins);
                    row = emptyNoiseRobustnessTable();
                    row(1,:) = {string(phantomName),string(noiseModel), ...
                        string(geometry),string(outputs(a).Algorithm), ...
                        cfg.noiseRobustnessSNRdB,measuredSNR,realization, ...
                        noiseSigma,metrics.PSNR,metrics.SSIM,metrics.Dice, ...
                        metrics.HD95MM,metrics.gCNR,metrics.ArtifactToTargetDB};
                    results = [results;row]; %#ok<AGROW>
                end
            end
        end
    end
end
end

%% ========================================================================
function tableOut = emptyNoiseRobustnessTable()
names = {'Phantom','NoiseModel','SensorGeometry','Algorithm', ...
    'RequestedSNR_dB','MeasuredSNR_dB','NoiseRealization','NominalNoiseSigma', ...
    'PSNR_ROI_dB','SSIM','Dice','HD95_mm','gCNR','ArtifactToTarget_dB'};
types = [{'string','string','string','string'},repmat({'double'},1,10)];
tableOut = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function [noisyData,sigma,measuredSNR] = addReviewerNoise(cleanData,snrDb,model,cfg)
signalRMS = rootMeanSquare(cleanData(:));
sigma = signalRMS*10^(-snrDb/20);
[M,N] = size(cleanData);
model = upper(string(model));
switch model
    case "IID_GAUSSIAN"
        noise = randn(M,N);
    case "CHANNEL_CORRELATED_GAUSSIAN"
        rho = min(max(cfg.channelNoiseCorrelation,0),0.95);
        common = randn(1,N);
        independent = randn(M,N);
        noise = sqrt(rho)*repmat(common,M,1)+sqrt(1-rho)*independent;
    case "HETEROSCEDASTIC_GAUSSIAN"
        scales = 1 + cfg.channelNoiseScaleSpread*linspace(-1,1,M).';
        noise = scales.*randn(M,N);
    case "LAPLACE_IID"
        u = rand(M,N)-0.5;
        noise = -sign(u).*log(max(1-2*abs(u),eps))/sqrt(2);
    otherwise
        error('Unknown reviewer noise model: %s',char(model));
end
noise = noise/rootMeanSquare(noise(:))*sigma;
noisyData = cleanData+noise;
measuredSNR = 20*log10(signalRMS/max(rootMeanSquare(noise(:)),eps));
end

%% ========================================================================
function savePSFProfileFigure(psf,rec,xVec,yVec,pointIndex,pointMM, ...
    geometry,algorithm,snrDb,outputDir)
fig = figure('Visible','off','Position',[100,100,1050,430]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;
imagesc(yVec*1e3,xVec*1e3,rec);
axis image; set(gca,'YDir','normal'); colorbar;
hold on; plot(pointMM(2),pointMM(1),'wx','MarkerSize',10,'LineWidth',1.5);
xlabel('y [mm]'); ylabel('x [mm]');
title('Normalized point reconstruction');
nexttile;
plot(psf.DistanceMM,psf.RadialProfile,'LineWidth',1.4);
hold on;
plot(psf.DistanceMM,psf.TangentialProfile,'--','LineWidth',1.4);
yline(0.5,':');
grid on;
xlabel('Distance from nominal absorber [mm]');
ylabel('Normalized amplitude');
legend('Radial','Tangential','Half maximum','Location','best');
title(sprintf(['FWHM_r %.3f mm, FWHM_t %.3f mm\n' ...
    'PSL_r %.1f dB, PSL_t %.1f dB'],psf.FWHMRadialMM, ...
    psf.FWHMTangentialMM,psf.PSLRadialDB,psf.PSLTangentialDB));
sgtitle(sprintf('Point %d (%.1f, %.1f mm) | %s | %s | SNR %s dB', ...
    pointIndex,pointMM(1),pointMM(2),char(geometry),char(algorithm), ...
    formatSNR(snrDb)));
fileName = sprintf('psf_p%02d_%s_%s_snr_%s.png',pointIndex, ...
    safeName(geometry),safeName(algorithm),safeName(formatSNR(snrDb)));
ensureDirectory(outputDir);
destination = fullfile(outputDir,fileName);
cleanupFigure = onCleanup(@() closeFigureIfValidV91(fig)); %#ok<NASGU>
try
    exportgraphics(fig,destination,'Resolution',180);
catch exportError
    ensureDirectory(fileparts(destination));
    try
        saveas(fig,destination);
    catch fallbackError
        warning('PACT:FigureExportFailed', ...
            'Could not save %s. exportgraphics: %s | saveas: %s', ...
            destination, exportError.message, fallbackError.message);
    end
end
end

%% ========================================================================
function [results,slopes,equivalence] = runMVComplexityBenchmark(cfg)
outputDir = fullfile(cfg.outputDir,'complexity_benchmark');
ensureDirectory(outputDir);
results = emptyComplexityTable();
equivalence = emptySyntheticEquivalenceTable();

rng(cfg.masterSeed+900009,'twister');
for M = cfg.complexitySensorCounts
    for L = cfg.complexitySubarrayLengths
        if L>floor(M/2)
            continue;
        end
        numPixels = cfg.complexityNumPixels;
        commonSignal = randn(1,numPixels);
        gain = 0.8+0.4*rand(M,1);
        delayed = gain*commonSignal + 0.6*randn(M,numPixels);

        localCfg = cfg;
        localCfg.mvSubarrayLength = L;
        localCfg.mvTilePixels = min(cfg.mvTilePixels,numPixels);
        context.isCircularArray = false;
        context.pixelMask = true(numPixels,1);
        context.noiseVariance = 0.36;
        context.numSensors = M;

        for w = 1:cfg.complexityWarmups
            spatiallySmoothedMVTiled(delayed,localCfg,context,"DIRECT");
            spatiallySmoothedMVTiled(delayed,localCfg,context,"SHERMAN_MORRISON");
        end

        directTimes = zeros(cfg.complexityRepetitions,1);
        smTimes = zeros(cfg.complexityRepetitions,1);
        directOutput = [];
        smOutput = [];
        for r = 1:cfg.complexityRepetitions
            tic;
            directOutput = spatiallySmoothedMVTiled(delayed,localCfg,context,"DIRECT");
            directTimes(r) = toc;
            tic;
            smOutput = spatiallySmoothedMVTiled(delayed,localCfg,context, ...
                "SHERMAN_MORRISON");
            smTimes(r) = toc;
        end

        relativeError = norm(directOutput-smOutput)/max(norm(directOutput),eps);
        maxDifference = max(abs(directOutput-smOutput));
        residual = meanInverseResidual(delayed,L,context.noiseVariance, ...
            cfg.mvDiagonalLoading,cfg.mvLoadingFloorRelative, ...
            cfg.inverseResidualSamples);

        K = M-L+1;
        memoryZMB = 8*L*K*min(cfg.mvTilePixels,numPixels)/1024^2;
        memoryPMB = 8*L^2*min(cfg.mvTilePixels,numPixels)/1024^2;

        results = [results;makeComplexityRows(M,L,K,numPixels,directTimes, ...
            smTimes,memoryZMB,memoryPMB)]; %#ok<AGROW>
        eqRow = emptySyntheticEquivalenceTable();
        eqRow(1,:) = {M,L,K,numPixels,relativeError,maxDifference,residual, ...
            median(directTimes),median(smTimes), ...
            median(directTimes)/max(median(smTimes),eps)};
        equivalence = [equivalence;eqRow]; %#ok<AGROW>

        fprintf(['Complexity: M=%3d, L=%2d | direct %.4f s | SM %.4f s | ' ...
            'rel. error %.3e\n'],M,L,median(directTimes),median(smTimes), ...
            relativeError);
    end
end
slopes = fitComplexitySlopes(results);
end

%% ========================================================================
function tableOut = emptyComplexityTable()
names = {'Method','NumSensors_M','SubarrayLength_L','NumSnapshots_K', ...
    'NumPixels','MedianRuntime_Sec','RuntimeIQR_Sec','MeanRuntime_Sec', ...
    'StdRuntime_Sec','EstimatedZMemory_MB','EstimatedStateMemory_MB'};
types = {'string','double','double','double','double','double','double', ...
    'double','double','double','double'};
tableOut = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function rows = makeComplexityRows(M,L,K,numPixels,directTimes,smTimes,zMB,pMB)
rows = emptyComplexityTable();
rows(1,:) = {"SS-MV Direct",M,L,K,numPixels,median(directTimes), ...
    percentileLinear(directTimes,75)-percentileLinear(directTimes,25), ...
    mean(directTimes),std(directTimes),zMB,zMB+8*L^2/1024^2};
rows(2,:) = {"SM-MV",M,L,K,numPixels,median(smTimes), ...
    percentileLinear(smTimes,75)-percentileLinear(smTimes,25), ...
    mean(smTimes),std(smTimes),zMB,zMB+pMB};
end

%% ========================================================================
function tableOut = emptySyntheticEquivalenceTable()
names = {'NumSensors_M','SubarrayLength_L','NumSnapshots_K','NumPixels', ...
    'RelativeL2Error','MaximumAbsoluteDifference','MeanInverseResidual', ...
    'DirectMedianRuntime_Sec','SMMedianRuntime_Sec','Speedup_DirectOverSM'};
types = repmat({'double'},1,numel(names));
tableOut = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);
end

%% ========================================================================
function residual = meanInverseResidual(delayed,L,noiseVariance,loadingFactor, ...
    floorRelative,numSamples)
[M,numPixels] = size(delayed);
K = M-L+1;
indices = unique(round(linspace(1,numPixels,min(numSamples,numPixels))));
values = zeros(numel(indices),1);
globalPower = mean(abs(delayed(:)).^2);
for i = 1:numel(indices)
    Z = buildSubarrayTensor(delayed,indices(i),L,K);
    Z = Z(:,:,1);
    R0 = (Z*Z')/K;
    loading = loadingFactor*(trace(R0)/L+noiseVariance)+ ...
        floorRelative*max(globalPower,eps);
    R = R0+loading*eye(L);
    P = eye(L)/loading;
    for k = 1:K
        u = Z(:,k)/sqrt(K);
        Pu = P*u;
        P = P-(Pu*Pu')/(1+u'*Pu);
    end
    values(i) = norm(R*P-eye(L),'fro')/sqrt(L);
end
residual = mean(values);
end

%% ========================================================================
function slopes = fitComplexitySlopes(results)
% Build the output table directly.  This avoids version-dependent MATLAB
% errors caused by assigning scalar structures containing mixed string and
% numeric fields into an initially empty structure array.
names = {'Method','SubarrayLength_L','LogLogSlope','Intercept', ...
    'RSquared','NumPoints'};
types = {'string','double','double','double','double','double'};
slopes = table('Size',[0,numel(names)],'VariableTypes',types, ...
    'VariableNames',names);

if isempty(results)
    return;
end

methods = unique(results.Method,'stable');
lengths = unique(results.SubarrayLength_L,'stable');
row = 0;
for m = 1:numel(methods)
    for l = 1:numel(lengths)
        idx = results.Method==methods(m) & ...
            results.SubarrayLength_L==lengths(l);
        if nnz(idx)<2
            continue;
        end

        sensorCounts = double(results.NumSensors_M(idx));
        runtimes = double(results.MedianRuntime_Sec(idx));
        valid = isfinite(sensorCounts) & sensorCounts>0 & ...
            isfinite(runtimes) & runtimes>0;
        if nnz(valid)<2
            continue;
        end

        x = log(sensorCounts(valid));
        y = log(runtimes(valid));
        coefficients = polyfit(x,y,1);
        fitted = polyval(coefficients,x);
        ssResidual = sum((y-fitted).^2);
        ssTotal = sum((y-mean(y)).^2);
        if ssTotal<=eps
            rSquared = NaN;
        else
            rSquared = 1-ssResidual/ssTotal;
        end

        row = row+1;
        slopes(row,:) = {string(methods(m)),double(lengths(l)), ...
            coefficients(1),coefficients(2),rSquared,nnz(valid)};
    end
end
end

%% ========================================================================
function tf = shouldSaveReconstruction(snrDb,realization,cfg)
tf = realization==cfg.saveRealization && anySNRMatch(snrDb,cfg.saveSNRdB);
end

%% ========================================================================
function saveReconstructionArray(recNorm,rawImage,phantom,geometry,algorithm, ...
    snrDb,realization,imageDir,cfg)
if ~cfg.saveRawReconstructions
    return;
end
fileName = sprintf('recon_%s_%s_snr_%s_rep_%02d.mat', ...
    safeName(geometry),safeName(algorithm),safeName(formatSNR(snrDb)),realization);
save(fullfile(imageDir,fileName),'recNorm','rawImage','phantom','geometry', ...
    'algorithm','snrDb','realization','-v7.3');
end

%% ========================================================================
function saveCaseMontage(gt,entries,xVec,yVec,phantom,geometry,snrDb, ...
    realization,imageDir)
numEntries = numel(entries);
numColumns = 4;
numRows = ceil((numEntries+1)/numColumns);
fig = figure('Visible','off','Position',[50,50,360*numColumns,320*numRows]);
t = tiledlayout(numRows,numColumns,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(yVec*1e3,xVec*1e3,gt);
axis image; set(gca,'YDir','normal'); colorbar;
title('Ground truth'); xlabel('y [mm]'); ylabel('x [mm]');

for k = 1:numEntries
    nexttile;
    imagesc(yVec*1e3,xVec*1e3,entries(k).Image);
    axis image; set(gca,'YDir','normal'); colorbar;
    title(sprintf('%s\nPSNR %.2f, gCNR %.3f', ...
        char(entries(k).Algorithm),entries(k).Metrics.PSNR, ...
        entries(k).Metrics.gCNR));
    xlabel('y [mm]'); ylabel('x [mm]');
end

title(t,sprintf('%s | %s | SNR %s dB | realization %d', ...
    char(phantom),char(geometry),formatSNR(snrDb),realization), ...
    'FontWeight','bold');
fileName = sprintf('montage_%s_snr_%s_rep_%02d.png',safeName(geometry), ...
    safeName(formatSNR(snrDb)),realization);
exportgraphics(fig,fullfile(imageDir,fileName),'Resolution',180);
close(fig);
end

%% ========================================================================
function text = safeName(value)
text = lower(char(string(value)));
text = regexprep(text,'[^a-zA-Z0-9]+','_');
text = regexprep(text,'^_+|_+$','');
if isempty(text)
    text = 'value';
end
end

%% ========================================================================
function saveSummaryFigures(summary,outputDir)
if isempty(summary)
    return;
end
if any(summary.EvaluationCohort=="ADAPTIVE_COMMON")
    summary = summary(summary.EvaluationCohort=="ADAPTIVE_COMMON",:);
end
figureDir = fullfile(outputDir,'paper_figures');
ensureDirectory(figureDir);
metrics = { ...
    'PSNR_ROI_dB_Mean','PSNR_ROI_dB_Std','PSNR [dB]','psnr_vs_snr.png'; ...
    'SSIM_Mean','SSIM_Std','SSIM','ssim_vs_snr.png'; ...
    'gCNR_Mean','gCNR_Std','gCNR','gcnr_vs_snr.png'; ...
    'HD95_mm_Mean','HD95_mm_Std','HD95 [mm]','hd95_vs_snr.png'; ...
    'ArtifactToTarget_dB_Mean','ArtifactToTarget_dB_Std', ...
        'Artifact-to-target power [dB]','artifact_power_vs_snr.png'};

for m = 1:size(metrics,1)
    fig = figure('Visible','off','Position',[100,100,1100,460]);
    tiledlayout(1,numel(unique(summary.SensorGeometry)),'Padding','compact');
    geometries = unique(summary.SensorGeometry,'stable');
    for g = 1:numel(geometries)
        nexttile;
        hold on; grid on;
        block = summary(summary.SensorGeometry==geometries(g),:);
        algorithms = unique(block.Algorithm,'stable');
        for a = 1:numel(algorithms)
            rows = block.Algorithm==algorithms(a);
            snr = block.RequestedSNR_dB(rows);
            [x,order] = plotSNRCoordinates(snr);
            y = block.(metrics{m,1})(rows);
            errorValues = block.(metrics{m,2})(rows);
            errorbar(x(order),y(order),errorValues(order),'-o', ...
                'DisplayName',char(algorithms(a)),'LineWidth',1.0);
        end
        [tickCoordinates,tickLabels] = snrAxisTicks(block.RequestedSNR_dB);
        xticks(tickCoordinates);
        xticklabels(tickLabels);
        xlabel('Requested SNR [dB]'); ylabel(metrics{m,3});
        title(strrep(char(geometries(g)),'_','\_'));
        legend('Location','bestoutside');
    end
    exportgraphics(fig,fullfile(figureDir,metrics{m,4}),'Resolution',180);
    close(fig);
end

% Runtime bar chart aggregated over all conditions.
algorithms = unique(summary.Algorithm,'stable');
medianRuntime = zeros(numel(algorithms),1);
for a = 1:numel(algorithms)
    medianRuntime(a) = median(summary.TotalRuntime_Sec_Median( ...
        summary.Algorithm==algorithms(a)),'omitnan');
end
fig = figure('Visible','off','Position',[100,100,1100,480]);
bar(categorical(algorithms,algorithms),medianRuntime);
ylabel('Median total runtime [s]');
grid on;
title('Pipeline runtime by reconstruction method');
xtickangle(35);
exportgraphics(fig,fullfile(figureDir,'runtime_by_algorithm.png'),'Resolution',180);
close(fig);
end

%% ========================================================================
function [coordinates,order] = plotSNRCoordinates(snr)
snr = double(snr(:));
finiteValues = snr(isfinite(snr));
if isempty(finiteValues)
    replacement = 0;
else
    replacement = max(finiteValues)+5;
end
coordinates = snr;
coordinates(isinf(coordinates)) = replacement;
[~,order] = sort(coordinates);
end

%% ========================================================================

function [coordinates,labels] = snrAxisTicks(snr)
values = unique(double(snr(:)),'stable');
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    replacement = 0;
else
    replacement = max(finiteValues)+5;
end
coordinates = values;
coordinates(isinf(coordinates)) = replacement;
[coordinates,order] = sort(coordinates);
values = values(order);
labels = strings(size(values));
for k = 1:numel(values)
    if isinf(values(k))
        labels(k) = "Inf";
    else
        labels(k) = string(values(k));
    end
end
end

%% ========================================================================
function saveComplexityFigures(results,slopes,outputDir)
figureDir = fullfile(outputDir,'complexity_benchmark');
ensureDirectory(figureDir);
lengths = unique(results.SubarrayLength_L,'stable');
methods = unique(results.Method,'stable');

for l = 1:numel(lengths)
    fig = figure('Visible','off','Position',[100,100,760,520]);
    hold on; grid on;
    for m = 1:numel(methods)
        idx = results.SubarrayLength_L==lengths(l) & results.Method==methods(m);
        block = sortrows(results(idx,:),'NumSensors_M');
        loglog(block.NumSensors_M,block.MedianRuntime_Sec,'-o', ...
            'LineWidth',1.4,'DisplayName',char(methods(m)));
    end
    xlabel('Number of receive elements M');
    ylabel('Median runtime [s]');
    title(sprintf('MV complexity benchmark, L=%d',lengths(l)));
    legend('Location','best');
    fileName = sprintf('mv_complexity_L_%02d.png',lengths(l));
    exportgraphics(fig,fullfile(figureDir,fileName),'Resolution',180);
    close(fig);
end

if ~isempty(slopes)
    writetable(slopes,fullfile(figureDir,'complexity_slopes_copy.csv'));
end
end

%% ========================================================================
function writeOutputDictionary(outputDir)
fileName = fullfile(outputDir,'OUTPUT_DICTIONARY.txt');
fid = fopen(fileName,'w');
if fid<0
    warning('Could not write output dictionary.');
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'PACT v9.0 output dictionary\n\n');
fprintf(fid,'all_reconstruction_results.csv\n');
fprintf(fid,'  One row per phantom, geometry, SNR, realization, and algorithm.\n\n');
fprintf(fid,'statistical_summary.csv\n');
fprintf(fid,'  Phantom-cluster mean, SD, median, IQR, and approximate 95%% CI by condition.\n\n');
fprintf(fid,'paired_comparisons_vs_DAS.csv\n');
fprintf(fid,'  Paired bootstrap intervals, sign-flip permutation p-values, and FDR q-values.\n\n');
fprintf(fid,'algorithm_rankings.csv\n');
fprintf(fid,'  Rank aggregation from normalized structural metrics only; not a general superiority claim.\n\n');
fprintf(fid,'metric_leaders_by_condition.csv\n');
fprintf(fid,'  Metric-specific extrema under the declared normalized evaluation, not an overall winner.\n\n');
fprintf(fid,'mv_image_equivalence.csv\n');
fprintf(fid,'  Image-domain SS-MV versus SM-MV numerical agreement.\n\n');
fprintf(fid,'mv_synthetic_equivalence.csv\n');
fprintf(fid,'  Controlled MV output error and inverse residual.\n\n');
fprintf(fid,'mv_complexity_results.csv and mv_complexity_slopes.csv\n');
fprintf(fid,'  Repeated runtimes, memory estimates, and empirical log-log slopes.\n\n');
fprintf(fid,'resolution_results.csv\n');
fprintf(fid,'  Radial/tangential FWHM, PSL, localization error, and runtime.\n\n');
fprintf(fid,'parameter_sensitivity.csv\n');
fprintf(fid,'  Multi-phantom, multi-geometry, multi-SNR parameter sensitivity.\n\n');
fprintf(fid,'delay_mismatch_results.csv\n');
fprintf(fid,'  Reconstruction-delay robustness under homogeneous and perturbed sound-speed maps.\n\n');
fprintf(fid,'noise_robustness_results.csv\n');
fprintf(fid,'  -5 dB robustness under IID, correlated, heteroscedastic, and Laplace noise.\n\n');
fprintf(fid,'Cohorts: CHEAP_ALL_PHANTOMS and ADAPTIVE_COMMON prevent unequal-cohort rankings.\n');
fprintf(fid,'Inferential tests use phantom-clustered differences, not individual noise repeats.\n');
fprintf(fid,'Important: no file in this benchmark reports iterative convergence.\n');
end

%% ========================================================================
function [karray,sensorX,sensorY,limitedIndices] = ...
    makeFiniteApertureArrayV90(kgrid,cfg)
theta = 2*pi*(0:cfg.numSensorsFull-1)/cfg.numSensorsFull;
sensorX = cfg.sensorRadius*cos(theta);
sensorY = cfg.sensorRadius*sin(theta);

karray = kWaveArray('BLITolerance',cfg.arrayBLITolerance, ...
    'UpsamplingRate',cfg.arrayUpsamplingRate);
for s = 1:cfg.numSensorsFull
    tangent = [-sin(theta(s)),cos(theta(s))];
    center = [sensorX(s),sensorY(s)];
    startPoint = center-0.5*cfg.elementWidth*tangent;
    endPoint = center+0.5*cfg.elementWidth*tangent;
    karray.addLineElement(startPoint,endPoint);
end
mask = karray.getArrayBinaryMask(kgrid);
if ~any(mask(:))
    error('kWaveArray generated an empty sensor mask.');
end

% A contiguous half-ring centered on the positive y direction.
centerIndex = round(cfg.numSensorsFull/4)+1;
halfWidth = floor(cfg.numSensorsLimited/2);
raw = (centerIndex-halfWidth):(centerIndex+cfg.numSensorsLimited-halfWidth-1);
limitedIndices = mod(raw-1,cfg.numSensorsFull)+1;
end

%% ========================================================================
function [phantoms,names] = generateArticlePhantoms(X,Y)
[native,nativeNames] = generateBasePhantomsV90(X,Y);
vascular = makeVascularTreePhantom(X,Y);
phantoms = [native;{vascular}];
names = [nativeNames,"Vascular_Tree"];
end

%% ========================================================================
function [phantoms,names] = generateBasePhantomsV90(X,Y)
names = ["Circle","Ring","Two_Circles","Ellipse","Rectangle", ...
    "Triangle","Shepp_Logan","Complex","Medical","Letter_A"];
phantoms = cell(numel(names),1);

phantoms{1} = finalizeContrast(double(hypot(X,Y)<=2.8e-3));
r = hypot(X,Y);
phantoms{2} = finalizeContrast(double(r>=2.3e-3 & r<=3.4e-3));
p = 0.75*double(hypot(X+2.5e-3,Y)<=1.8e-3)+ ...
    1.00*double(hypot(X-2.4e-3,Y-0.7e-3)<=1.2e-3);
phantoms{3} = finalizeContrast(p);
phantoms{4} = finalizeContrast(double((X/4.0e-3).^2+(Y/2.2e-3).^2<=1));
phantoms{5} = finalizeContrast(double(abs(X)<=3.8e-3 & abs(Y)<=2.2e-3));
vertexX = [-4.2e-3,4.2e-3,0];
vertexY = [3.5e-3,3.5e-3,-4.5e-3];
phantoms{6} = finalizeContrast(double(inpolygon(X,Y,vertexX,vertexY)));
phantoms{7} = finalizeContrast(customSheppLogan(X,Y));

p = 0.05*double(hypot(X,Y)<=8.7e-3);
p = max(p,0.55*double(hypot(X,Y)<=5.2e-3));
p = max(p,0.78*double(abs(X)<=2.8e-3 & abs(Y)<=2.8e-3));
p = max(p,1.00*double(hypot(X-5.0e-3,Y+3.5e-3)<=1.2e-3));
p = max(p,0.68*double(abs(X+5.7e-3)<=0.45e-3 & abs(Y-1.5e-3)<=3.5e-3));
tri = inpolygon(X,Y,[-5.4e-3,-2.4e-3,-3.9e-3],[4.8e-3,4.8e-3,1.8e-3]);
p = max(p,0.88*double(tri));
phantoms{8} = finalizeContrast(gaussianSmooth2D(p,0.7));

support = hypot(X,Y)<=8.7e-3;
texture = 0.04+0.012*sin(2*pi*X/3.2e-3).*cos(2*pi*Y/4.0e-3);
p = max(texture,0).*support;
tumorCenterX = -1.5e-3;
tumorCenterY = -0.8e-3;
tumorRadius = 3.0e-3;
tumorDistance = hypot(X-tumorCenterX,Y-tumorCenterY);
tumor = tumorDistance<=tumorRadius;
tumorIntensity = 0.65+0.35*(1-tumorDistance/tumorRadius);
p(tumor) = max(p(tumor),tumorIntensity(tumor));
necrosis = hypot(X-(tumorCenterX+0.5e-3),Y-(tumorCenterY+0.4e-3))<=1.25e-3;
p(necrosis) = 0.10;
v1 = distanceToSegment(X,Y,[3.8e-3,-4.8e-3],[3.8e-3,4.8e-3])<=0.22e-3;
v2 = distanceToSegment(X,Y,[-4.8e-3,-3.2e-3],[4.8e-3,-3.2e-3])<=0.22e-3;
v3 = distanceToSegment(X,Y,[3.8e-3,-3.2e-3],[6.5e-3,-1.8e-3])<=0.18e-3;
p(v1) = max(p(v1),0.95);
p(v2) = max(p(v2),0.85);
p(v3) = max(p(v3),0.75);
calcCenters = [-5.8e-3,4.6e-3;5.6e-3,-2.5e-3;-3.6e-3,-5.5e-3;2.8e-3,5.6e-3];
for k = 1:size(calcCenters,1)
    calc = hypot(X-calcCenters(k,1),Y-calcCenters(k,2))<=0.40e-3;
    p(calc) = 1;
end
p(~support) = 0;
phantoms{9} = finalizeContrast(gaussianSmooth2D(p,0.55));

leftLeg = distanceToSegment(X,Y,[4.8e-3,0],[-4.5e-3,-4.0e-3])<=0.45e-3;
rightLeg = distanceToSegment(X,Y,[4.8e-3,0],[-4.5e-3,4.0e-3])<=0.45e-3;
crossbar = distanceToSegment(X,Y,[0.2e-3,-2.1e-3],[0.2e-3,2.1e-3])<=0.38e-3;
phantoms{10} = finalizeContrast(gaussianSmooth2D(double(leftLeg|rightLeg|crossbar),0.45));
end

%% ========================================================================
function p = makeVascularTreePhantom(X,Y)
p = zeros(size(X));
segments = [ ...
    -6.8,-0.2, 5.8, 0.2, 0.34,1.00; ...
    -2.8, 0.0, 1.5, 4.8, 0.25,0.90; ...
    -2.2, 0.0, 0.4,-4.9, 0.24,0.88; ...
     0.8, 0.1, 5.2, 3.7, 0.22,0.85; ...
     1.2, 0.1, 5.7,-3.5, 0.20,0.82; ...
    -0.3, 2.7, 2.0, 6.1, 0.16,0.78; ...
    -0.4, 2.6,-3.2, 5.7, 0.15,0.76; ...
    -0.7,-2.7, 2.4,-6.0, 0.15,0.74; ...
    -0.8,-2.6,-3.7,-5.5, 0.14,0.72; ...
     3.6, 2.4, 6.7, 4.9, 0.12,0.68; ...
     3.8,-2.0, 6.8,-4.6, 0.11,0.66; ...
    -2.5, 4.8,-5.4, 6.7, 0.10,0.62; ...
    -2.8,-4.5,-5.8,-6.6, 0.10,0.60];
for k = 1:size(segments,1)
    p1 = segments(k,1:2)*1e-3;
    p2 = segments(k,3:4)*1e-3;
    radius = segments(k,5)*1e-3;
    intensity = segments(k,6);
    mask = distanceToSegment(X,Y,p1,p2)<=radius;
    p(mask) = max(p(mask),intensity);
end
% Add a few small vessel cross-sections.
centers = [-4.2,2.0;2.8,5.2;5.7,0.8;-2.4,-6.2]*1e-3;
for k = 1:size(centers,1)
    mask = hypot(X-centers(k,1),Y-centers(k,2))<=0.22e-3;
    p(mask) = max(p(mask),0.80);
end
p(hypot(X,Y)>8.4e-3) = 0;
p = finalizeContrast(gaussianSmooth2D(p,0.35));
end

%% ========================================================================
function background = buildClinicalBackground(X, Y, cfg)
% Fixed illustrative anatomy used for every optical target. This avoids
% changing the acoustic inverse problem when only the optical target changes.
r = hypot(X, Y);
bodyMask = r <= cfg.bodyRadius;
fatMask = bodyMask & r >= 7.3e-3;
innerMask = bodyMask & ~fatMask;

% Fixed low-frequency texture, independent of phantom index.
rng(88031, 'twister');
texture = gaussianSmooth2D(randn(size(X)), 9);
texture = texture / (max(abs(texture(:))) + eps);

% A few denser connective/glandular regions.
connective1 = ((X+2.1e-3)/2.6e-3).^2 + ((Y-1.2e-3)/1.5e-3).^2 <= 1;
connective2 = ((X-2.8e-3)/1.8e-3).^2 + ((Y+2.0e-3)/2.7e-3).^2 <= 1;
connectiveMask = innerMask & (connective1 | connective2);

% Acoustic maps. Outside the body is water coupling medium.
soundSpeed = 1480 * ones(size(X));                    % [m/s]
density = 998 * ones(size(X));                        % [kg/m^3]
alphaCoeff = 0.002 * ones(size(X));                   % dB/(MHz^y cm)

soundSpeed(fatMask) = 1450 + 8*texture(fatMask);
density(fatMask) = 930 + 8*texture(fatMask);
alphaCoeff(fatMask) = 0.55 + 0.05*texture(fatMask);

soundSpeed(innerMask) = 1540 + 18*texture(innerMask);
density(innerMask) = 1000 + 12*texture(innerMask);
alphaCoeff(innerMask) = 0.75 + 0.08*texture(innerMask);

soundSpeed(connectiveMask) = 1580 + 12*texture(connectiveMask);
density(connectiveMask) = 1040 + 10*texture(connectiveMask);
alphaCoeff(connectiveMask) = 0.95 + 0.08*texture(connectiveMask);

% Enforce physically meaningful ranges.
soundSpeed = min(max(soundSpeed, 1430), 1610);
density = min(max(density, 900), 1080);
alphaCoeff = min(max(alphaCoeff, 0), 1.2);

% Baseline optical maps used by the fluence solver.
muABase = zeros(size(X));                              % [1/m]
muSpBase = zeros(size(X));                             % reduced scattering [1/m]
gammaBase = zeros(size(X));                            % dimensionless

muABase(fatMask) = 5.0;
muSpBase(fatMask) = 800;
gammaBase(fatMask) = 0.11;

muABase(innerMask) = 8.0;
muSpBase(innerMask) = 1100;
gammaBase(innerMask) = 0.13;

muABase(connectiveMask) = 11.0;
muSpBase(connectiveMask) = 1300;
gammaBase(connectiveMask) = 0.14;

background.bodyMask = bodyMask;
background.fatMask = fatMask;
background.connectiveMask = connectiveMask;
background.soundSpeed = single(soundSpeed);
background.density = single(density);
background.alphaCoeff = single(alphaCoeff);
background.muABase = muABase;
background.muSpBase = muSpBase;
background.gammaBase = gammaBase;
end

%% ========================================================================

%% ========================================================================
function optical = buildOpticalAndPressureMaps(contrast, background, X, Y, cfg)
contrast = min(max(double(contrast), 0), 1);
contrast(~background.bodyMask) = 0;

muA = background.muABase + cfg.targetMuAContrast * contrast;
muSp = background.muSpBase .* (1 + cfg.targetMuSpFraction * contrast);
gamma = background.gammaBase .* (1 + cfg.targetGruneisenFraction * contrast);

% Avoid zeros inside the tissue and keep outside coupling water excluded
% from the optical solve.
muA(background.bodyMask) = max(muA(background.bodyMask), 0.1);
muSp(background.bodyMask) = max(muSp(background.bodyMask), 100);
gamma(background.bodyMask) = max(gamma(background.bodyMask), 0.01);

switch upper(char(cfg.fluenceModel))
    case 'DIFFUSION'
        fluence = solveDiffuseFluence2D(muA, muSp, background.bodyMask, ...
            X, Y, cfg);
    otherwise
        error('Unsupported fluenceModel: %s', char(cfg.fluenceModel));
end

% Dimensional thermoelastic initial pressure.
% mu_a [1/m] * fluence [J/m^2] = absorbed energy density [J/m^3 = Pa].
p0Pa = gamma .* muA .* fluence;
p0Pa(~background.bodyMask) = 0;

optical.contrast = contrast;
optical.muA_mInv = muA;
optical.muSp_mInv = muSp;
optical.gruneisen = gamma;
optical.fluenceJm2 = fluence;
optical.absorbedEnergyJm3 = muA .* fluence;
optical.p0Pa = p0Pa;
end

%% ========================================================================

%% ========================================================================
function fluence = solveDiffuseFluence2D(muA, muSp, tissueMask, X, Y, cfg)
% Solve the steady diffusion equation
%   -div(D grad Phi) + mu_a Phi = 0
% using prescribed illumination on the tissue boundary.
%
% This is a first-order approximation for highly scattering tissue. It is
% not a replacement for radiative transport or Monte Carlo near sources,
% boundaries, or low-scattering regions.

[Nx, Ny] = size(muA);
D = zeros(Nx, Ny);
D(tissueMask) = 1 ./ (3 * (muA(tissueMask) + muSp(tissueMask))); % [m]

boundaryMask = bwperim(tissueMask, 8);
interiorMask = tissueMask & ~boundaryMask;

boundaryFluence = makeBoundaryIllumination(boundaryMask, X, Y, cfg);
fluence = zeros(Nx, Ny);
fluence(boundaryMask) = boundaryFluence(boundaryMask);

unknownLinear = find(interiorMask);
numUnknown = numel(unknownLinear);
if numUnknown == 0
    error('No interior tissue pixels are available for the fluence solve.');
end

unknownIndex = zeros(Nx, Ny, 'uint32');
unknownIndex(unknownLinear) = uint32(1:numUnknown);

% At most five non-zero entries per equation.
rows = zeros(5*numUnknown, 1);
cols = zeros(5*numUnknown, 1);
vals = zeros(5*numUnknown, 1);
b = zeros(numUnknown, 1);
entry = 0;

neighborOffsets = [-1,0; 1,0; 0,-1; 0,1];
hx2 = cfg.dx^2;
hy2 = cfg.dy^2;

for u = 1:numUnknown
    [i, j] = ind2sub([Nx, Ny], unknownLinear(u));
    diagonal = muA(i,j);

    for q = 1:4
        ni = i + neighborOffsets(q,1);
        nj = j + neighborOffsets(q,2);

        if ni < 1 || ni > Nx || nj < 1 || nj > Ny
            continue;
        end

        if q <= 2
            h2 = hx2;
        else
            h2 = hy2;
        end

        if tissueMask(ni,nj)
            faceD = harmonicMean(D(i,j), D(ni,nj));
            coefficient = faceD / h2;
            diagonal = diagonal + coefficient;

            if interiorMask(ni,nj)
                entry = entry + 1;
                rows(entry) = u;
                cols(entry) = double(unknownIndex(ni,nj));
                vals(entry) = -coefficient;
            else
                b(u) = b(u) + coefficient * fluence(ni,nj);
            end
        else
            % Zero-fluence exterior boundary. Boundary pixels are already
            % excluded from the unknown set, so this term is usually absent.
            faceD = D(i,j);
            coefficient = faceD / h2;
            diagonal = diagonal + coefficient;
        end
    end

    entry = entry + 1;
    rows(entry) = u;
    cols(entry) = u;
    vals(entry) = diagonal;
end

A = sparse(rows(1:entry), cols(1:entry), vals(1:entry), ...
    numUnknown, numUnknown);
phiInterior = A \ b;

if any(~isfinite(phiInterior)) || any(phiInterior < -1e-9)
    error('The optical diffusion solve produced invalid fluence values.');
end

fluence(unknownLinear) = max(phiInterior, 0);
fluence(~tissueMask) = 0;
end

%% ========================================================================

%% ========================================================================
function boundaryFluence = makeBoundaryIllumination(boundaryMask, X, Y, cfg)
theta = atan2(Y, X);
profile = zeros(size(X));
sigma = deg2rad(cfg.illuminationAngularSigmaDeg);
sourceAngles = 2*pi*(0:cfg.numIlluminationFibers-1) / cfg.numIlluminationFibers;

for k = 1:numel(sourceAngles)
    delta = atan2(sin(theta-sourceAngles(k)), cos(theta-sourceAngles(k)));
    profile = profile + exp(-0.5*(delta/sigma).^2);
end

boundaryValues = profile(boundaryMask);
profile = profile / (max(boundaryValues) + eps);
profile = cfg.illuminationFloorFraction + ...
    (1-cfg.illuminationFloorFraction) * profile;

boundaryFluence = zeros(size(X));
boundaryFluence(boundaryMask) = cfg.surfaceFluenceJm2 * profile(boundaryMask);
end

%% ========================================================================

%% ========================================================================
function value = harmonicMean(a, b)
value = 2*a*b / (a+b+eps);
end

%% ========================================================================

%% ========================================================================
function sensorData = simulateKWaveFiniteArray(p0Pa, background, ...
    kgrid, karray, cfg)
medium.sound_speed = background.soundSpeed;
medium.sound_speed_ref = double(max(background.soundSpeed(:)));
medium.density = background.density;
medium.alpha_coeff = background.alphaCoeff;
medium.alpha_power = cfg.alphaPower;

source.p0 = single(p0Pa);

sensor.mask = karray.getArrayBinaryMask(kgrid);
sensor.frequency_response = [cfg.sensorCenterFrequency, ...
    cfg.sensorBandwidthPercent];

inputArgs = { ...
    'PMLInside', false, ...
    'PMLSize', 20, ...
    'PlotSim', false, ...
    'PlotPML', false, ...
    'DataCast', 'single', ...
    'DataRecast', true, ...
    'Smooth', [true, false, false]};

rawGridData = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputArgs{:});
sensorData = karray.combineSensorData(kgrid, rawGridData);

if size(sensorData,1) ~= cfg.numSensorsFull
    error('Expected %d element traces, but kWaveArray returned %d.', ...
        cfg.numSensorsFull, size(sensorData,1));
end
end

%% ========================================================================

%% ========================================================================
function travelTimes = computeTravelTimes(X, Y, xVec, yVec, ...
    sensorX, sensorY, soundSpeed, cfg)
numSensors = numel(sensorX);
numPixels = numel(X);
travelTimes = zeros(numSensors, numPixels, 'single');

switch upper(char(cfg.delayModel))
    case 'HOMOGENEOUS'
        for s = 1:numSensors
            travelTimes(s,:) = single(hypot(X(:)-sensorX(s), ...
                Y(:)-sensorY(s)) / cfg.homogeneousReconSoundSpeed);
        end

    case 'STRAIGHT_RAY_HETEROGENEOUS'
        pixelX = X(:).';
        pixelY = Y(:).';
        cFallback = cfg.homogeneousReconSoundSpeed;
        fractions = ((1:cfg.straightRaySamples)-0.5) / cfg.straightRaySamples;

        for s = 1:numSensors
            dxRay = pixelX - sensorX(s);
            dyRay = pixelY - sensorY(s);
            distance = hypot(dxRay, dyRay);
            meanSlowness = zeros(1, numPixels);

            for q = 1:numel(fractions)
                fraction = fractions(q);
                sampleX = sensorX(s) + fraction*dxRay;
                sampleY = sensorY(s) + fraction*dyRay;

                % Matrix rows correspond to xVec and columns to yVec.
                localC = interp2(yVec, xVec, double(soundSpeed), ...
                    sampleY, sampleX, 'linear', cFallback);
                meanSlowness = meanSlowness + 1 ./ localC;
            end

            meanSlowness = meanSlowness / numel(fractions);
            travelTimes(s,:) = single(distance .* meanSlowness);
        end

    otherwise
        error('Unknown delayModel: %s', char(cfg.delayModel));
end
end

%% ========================================================================

%% ========================================================================
function delayed = sampleDelayedData(sensorData, tArray, travelTimes)
[numSensors, ~] = size(sensorData);
numPixels = size(travelTimes, 2);
delayed = zeros(numSensors, numPixels, 'single');

for s = 1:numSensors
    values = interp1(tArray, sensorData(s,:), ...
        double(travelTimes(s,:)), 'linear', 0);
    delayed(s,:) = single(values);
end
end

%% ========================================================================

%% ========================================================================
function saveBackgroundFigure(background, xVec, yVec, outputDir)
fig = figure('Position', [100, 100, 1500, 760], 'Visible', 'off');

subplot(2,3,1);
imagesc(yVec*1e3, xVec*1e3, background.soundSpeed);
axis image; set(gca,'YDir','normal'); colorbar;
title('Sound speed [m/s]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,2);
imagesc(yVec*1e3, xVec*1e3, background.density);
axis image; set(gca,'YDir','normal'); colorbar;
title('Density [kg/m^3]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,3);
imagesc(yVec*1e3, xVec*1e3, background.alphaCoeff);
axis image; set(gca,'YDir','normal'); colorbar;
title('alpha_0 [dB/(MHz^y cm)]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,4);
imagesc(yVec*1e3, xVec*1e3, background.muABase);
axis image; set(gca,'YDir','normal'); colorbar;
title('Baseline \mu_a [1/m]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,5);
imagesc(yVec*1e3, xVec*1e3, background.muSpBase);
axis image; set(gca,'YDir','normal'); colorbar;
title('Baseline \mu_s'' [1/m]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,6);
imagesc(yVec*1e3, xVec*1e3, background.gammaBase);
axis image; set(gca,'YDir','normal'); colorbar;
title('Gruneisen parameter'); xlabel('y [mm]'); ylabel('x [mm]');

sgtitle('Fixed clinically motivated background maps');
exportgraphics(fig, fullfile(outputDir, 'background_physical_maps.png'), ...
    'Resolution', 180);
close(fig);
end

%% ========================================================================

%% ========================================================================
function savePhysicalMaps(optical, background, xVec, yVec, phantomDir, name)
save(fullfile(phantomDir, 'physical_maps.mat'), 'optical', 'background', '-v7.3');

fig = figure('Position', [100, 100, 1500, 760], 'Visible', 'off');

subplot(2,3,1);
imagesc(yVec*1e3, xVec*1e3, optical.muA_mInv);
axis image; set(gca,'YDir','normal'); colorbar;
title('\mu_a [1/m]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,2);
imagesc(yVec*1e3, xVec*1e3, optical.muSp_mInv);
axis image; set(gca,'YDir','normal'); colorbar;
title('\mu_s'' [1/m]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,3);
imagesc(yVec*1e3, xVec*1e3, optical.fluenceJm2);
axis image; set(gca,'YDir','normal'); colorbar;
title('Fluence [J/m^2]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,4);
imagesc(yVec*1e3, xVec*1e3, optical.gruneisen);
axis image; set(gca,'YDir','normal'); colorbar;
title('Gruneisen'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,5);
imagesc(yVec*1e3, xVec*1e3, optical.absorbedEnergyJm3);
axis image; set(gca,'YDir','normal'); colorbar;
title('Absorbed energy [J/m^3]'); xlabel('y [mm]'); ylabel('x [mm]');

subplot(2,3,6);
imagesc(yVec*1e3, xVec*1e3, optical.p0Pa);
axis image; set(gca,'YDir','normal'); colorbar;
title('Initial pressure p_0 [Pa]'); xlabel('y [mm]'); ylabel('x [mm]');

sgtitle(sprintf('%s: optical and photoacoustic maps', char(name)));
exportgraphics(fig, fullfile(phantomDir, 'optical_and_p0_maps.png'), ...
    'Resolution', 180);
close(fig);
end

%% ========================================================================

%% ========================================================================
function p = customSheppLogan(X, Y)
E = [ ...
     1.00, .6900, .9200,  0.000,  0.0000,   0;
    -0.80, .6624, .8740,  0.000, -0.0184,   0;
    -0.20, .1100, .3100,  0.220,  0.0000, -18;
    -0.20, .1600, .4100, -0.220,  0.0000,  18;
     0.10, .2100, .2500,  0.000,  0.3500,   0;
     0.10, .0460, .0460,  0.000,  0.1000,   0;
     0.10, .0460, .0460,  0.000, -0.1000,   0;
     0.10, .0460, .0230, -0.080, -0.6050,   0;
     0.10, .0230, .0230,  0.000, -0.6060,   0;
     0.10, .0230, .0460,  0.060, -0.6050,   0];

scale = 8.0e-3;
Xn = X / scale;
Yn = Y / scale;
p = zeros(size(X));
for k = 1:size(E,1)
    angle = deg2rad(E(k,6));
    xShift = Xn - E(k,4);
    yShift = Yn - E(k,5);
    xr =  xShift*cos(angle) + yShift*sin(angle);
    yr = -xShift*sin(angle) + yShift*cos(angle);
    mask = (xr/E(k,2)).^2 + (yr/E(k,3)).^2 <= 1;
    p(mask) = p(mask) + E(k,1);
end
p = max(p, 0);
end

%% ========================================================================

%% ========================================================================
function d = distanceToSegment(X, Y, p1, p2)
vx = p2(1) - p1(1);
vy = p2(2) - p1(2);
wx = X - p1(1);
wy = Y - p1(2);
t = (wx*vx + wy*vy) / (vx^2 + vy^2 + eps);
t = min(max(t, 0), 1);
projectionX = p1(1) + t*vx;
projectionY = p1(2) + t*vy;
d = hypot(X-projectionX, Y-projectionY);
end

%% ========================================================================

%% ========================================================================
function smoothed = gaussianSmooth2D(image, sigmaPixels)
if sigmaPixels <= 0
    smoothed = image;
    return;
end
radius = max(1, ceil(3*sigmaPixels));
x = -radius:radius;
g = exp(-(x.^2)/(2*sigmaPixels^2));
g = g / sum(g);
smoothed = conv2(conv2(image, g, 'same'), g.', 'same');
end

%% ========================================================================

%% ========================================================================
function p = finalizeContrast(p)
p = max(real(p), 0);
maximum = max(p(:));
if maximum <= eps
    error('Generated optical target is empty.');
end
p = p / maximum;
end
