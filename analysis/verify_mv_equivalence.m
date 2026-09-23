function S = verify_mv_equivalence(mainCsv,equivCsv,outCsv)
% Verify the stored full-ring SM-MV/DAS metric identity and the independent
% direct-SS-MV vs Sherman-Morrison solver equivalence. No reconstruction is run.

if nargin<1 || isempty(mainCsv)
    mainCsv = 'all_reconstruction_results_PARTIAL.csv';
end
if nargin<2 || isempty(equivCsv)
    equivCsv = 'mv_image_equivalence_PARTIAL.csv';
end
if nargin<3 || isempty(outCsv)
    outCsv = 'MV_fullring_collapse_validation.csv';
end

T = readtable(mainCsv,'TextType','string');
E = readtable(equivCsv,'TextType','string');

D = T(T.SensorGeometry=="CIRCULAR_FULL" & T.Algorithm=="DAS",:);
M = T(T.SensorGeometry=="CIRCULAR_FULL" & T.Algorithm=="SM-MV",:);

J = innerjoin(M,D,'Keys', ...
    {'Phantom','SensorGeometry','RequestedSNR_dB','NoiseRealization'}, ...
    'LeftVariables',{'Phantom','SensorGeometry','RequestedSNR_dB','NoiseRealization', ...
    'RMSE_ROI','PSNR_ROI_dB','SSIM','Dice','HD95_mm','gCNR','ArtifactToTarget_dB'}, ...
    'RightVariables',{'RMSE_ROI','PSNR_ROI_dB','SSIM','Dice','HD95_mm','gCNR','ArtifactToTarget_dB'});

metrics = ["RMSE_ROI","PSNR_ROI_dB","SSIM","Dice","HD95_mm","gCNR","ArtifactToTarget_dB"];
for k=1:numel(metrics)
    % innerjoin suffixes duplicate variables as _Tleft/_Tright in recent MATLAB
    names = string(J.Properties.VariableNames);
    candidates = names(startsWith(names,metrics(k)));
    if numel(candidates)~=2
        error('Could not resolve joined columns for %s.',metrics(k));
    end
    J.("Delta_"+metrics(k)) = J.(candidates(1)) - J.(candidates(2));
end

keep = ["Phantom","SensorGeometry","RequestedSNR_dB","NoiseRealization", ...
    "Delta_RMSE_ROI","Delta_PSNR_ROI_dB","Delta_SSIM","Delta_Dice", ...
    "Delta_HD95_mm","Delta_gCNR","Delta_ArtifactToTarget_dB"];
writetable(J(:,keep),outCsv);

fprintf('Full-ring stored SM-MV/DAS validation cases: %d\n',height(J));
fprintf('Max |Delta PSNR| = %.3g dB\n',max(abs(J.Delta_PSNR_ROI_dB)));
fprintf('Max |Delta SSIM| = %.3g\n',max(abs(J.Delta_SSIM)));
fprintf('Max |Delta Dice| = %.3g\n',max(abs(J.Delta_Dice)));
fprintf('Direct vs Sherman-Morrison cases: %d\n',height(E));
fprintf('Max relative L2 error = %.6g\n',max(E.RelativeL2Error));
fprintf('Max absolute difference = %.6g\n',max(E.MaximumAbsoluteDifference));
fprintf('Min image correlation = %.15g\n',min(E.ImageCorrelation));

S = struct();
S.fullRingPairs = J(:,keep);
S.directVsSM = E;
end
