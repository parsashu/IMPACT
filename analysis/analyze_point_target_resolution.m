function outputs = analyze_point_target_resolution(csvFile,outputDir)
% analyze_point_target_resolution
% Post-process the completed v9.7 point-target benchmark.
%
% Statistical convention:
%   five noise realizations are averaged within each fixed target location
%   first. Spatial summaries are then descriptive across the nine prescribed
%   locations; they are not treated as an inferential random sample.
%
% Outputs:
%   resolution_9point_point_means.csv
%   resolution_9point_spatial_summary.csv
%   Table_resolution_9point_20dB.csv
%
% The manuscript figures supplied with the revision are generated from these
% same point-level means.

if nargin<1 || isempty(csvFile)
    csvFile = fullfile(pwd,'RESOLUTION_9POINT_REVISION', ...
        'resolution_9point_results.csv');
end
if nargin<2 || isempty(outputDir)
    outputDir = fullfile(pwd,'RESOLUTION_9POINT_ANALYSIS_MATLAB');
end
if ~isfolder(outputDir), mkdir(outputDir); end

T = readtable(csvFile,'TextType','string');

assert(height(T)==1620, ...
    'Expected the completed 1620-row v9.7 table; found %d rows.',height(T));
assert(numel(unique(T.PointIndex))==9,'Expected nine target locations.');

metrics = ["FWHM_Radial_mm","FWHM_Tangential_mm", ...
    "PSL_Radial_dB","PSL_Tangential_dB","PeakLocalizationError_mm"];

% Average the five repeated noise realizations within each fixed position.
G = findgroups(T.PointIndex,T.PointX_mm,T.PointY_mm,T.SensorGeometry, ...
    T.Algorithm,T.RequestedSNR_dB);

P = table();
[P.PointIndex,P.PointX_mm,P.PointY_mm,P.SensorGeometry, ...
 P.Algorithm,P.RequestedSNR_dB] = splitapply(@firstTuple, ...
    T.PointIndex,T.PointX_mm,T.PointY_mm,T.SensorGeometry, ...
    T.Algorithm,T.RequestedSNR_dB,G);

for k=1:numel(metrics)
    P.(metrics(k)) = splitapply(@(x)mean(x,'omitnan'),T.(metrics(k)),G);
end
writetable(P,fullfile(outputDir,'resolution_9point_point_means.csv'));

% Descriptive spatial summary after realization averaging.
G2 = findgroups(P.SensorGeometry,P.Algorithm,P.RequestedSNR_dB);
S = table();
[S.SensorGeometry,S.Algorithm,S.RequestedSNR_dB] = ...
    splitapply(@firstTuple,P.SensorGeometry,P.Algorithm,P.RequestedSNR_dB,G2);
S.NumPositions = splitapply(@numel,P.PointIndex,G2);

for k=1:numel(metrics)
    x=P.(metrics(k));
    S.(metrics(k)+"_Mean") = splitapply(@(v)mean(v,'omitnan'),x,G2);
    S.(metrics(k)+"_SD")   = splitapply(@(v)std(v,0,'omitnan'),x,G2);
    S.(metrics(k)+"_Min")  = splitapply(@(v)min(v,[],'omitnan'),x,G2);
    S.(metrics(k)+"_Max")  = splitapply(@(v)max(v,[],'omitnan'),x,G2);
end
writetable(S,fullfile(outputDir,'resolution_9point_spatial_summary.csv'));

selectedFull = ["DAS","CF-DAS","SCF-DAS","DMAS","DASDSF"];
selectedLimited = ["DAS","SM-MV","CF-SM-MV","DASDSF-SM-MV","NCVDAS-SM-MV"];
K = S(S.RequestedSNR_dB==20 & ...
    ((S.SensorGeometry=="CIRCULAR_FULL" & ismember(S.Algorithm,selectedFull)) | ...
     (S.SensorGeometry=="LIMITED_VIEW_ARC" & ismember(S.Algorithm,selectedLimited))),:);
writetable(K,fullfile(outputDir,'Table_resolution_9point_20dB.csv'));

fprintf('Completed point-target rows: %d\n',height(T));
fprintf('Point-level condition means: %d\n',height(P));
fprintf('Spatial summary rows: %d\n',height(S));

outputs = struct('pointMeans',P,'spatialSummary',S,'table20dB',K, ...
    'outputDir',string(outputDir));
end

function varargout = firstTuple(varargin)
varargout = cell(size(varargin));
for k=1:nargin
    varargout{k}=varargin{k}(1);
end
end
