function PACT_collect_environment_info(outputFile)
% Record MATLAB, toolbox, and k-Wave provenance for the public release.
if nargin < 1 || isempty(outputFile)
    outputFile = 'ENVIRONMENT.txt';
end
fid = fopen(outputFile,'w');
if fid < 0, error('Cannot open %s',outputFile); end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid,'Generated: %s\n',char(datetime('now')));
fprintf(fid,'MATLAB version: %s\n',version);
fprintf(fid,'MATLAB release: %s\n',version('-release'));
fprintf(fid,'Architecture: %s\n',computer('arch'));
fprintf(fid,'Platform: %s\n\n',computer);

fprintf(fid,'Installed MATLAB products:\n');
V = ver;
for k=1:numel(V)
    fprintf(fid,'  %s | Version %s | Release %s\n', ...
        V(k).Name,V(k).Version,V(k).Release);
end

fprintf(fid,'\nk-Wave provenance:\n');
fprintf(fid,'  kspaceFirstOrder2D: %s\n',which('kspaceFirstOrder2D'));
fprintf(fid,'  kWaveArray: %s\n',which('kWaveArray'));
if exist('getkWavePath','file')==2
    try
        fprintf(fid,'  getkWavePath: %s\n',getkWavePath());
    catch ME
        fprintf(fid,'  getkWavePath error: %s\n',ME.message);
    end
end
if exist('getComputerInfo','file')==2
    try
        fprintf(fid,'\ngetComputerInfo output:\n');
        info = evalc('getComputerInfo');
        fprintf(fid,'%s\n',info);
    catch ME
        fprintf(fid,'  getComputerInfo error: %s\n',ME.message);
    end
end
end
