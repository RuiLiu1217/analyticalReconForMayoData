% analytical Reconstruction demo
% author: Rui Liu

% Read Projection data
arch = computer;

if strcmp(arch, 'GLNXA64') % LINUX
    ProjectionDataPath = '/home/liurui/Desktop/AAPM/AAPM_Data/L067/L067/full_DICOM-CT-PD';
    dictionaryFile =  '/home/liurui/Desktop/AAPM/AAPM_Data/GeneralInfo/GeneralInfo/DICOM-CT-PD-dict_v8.txt';
elseif  strcmp(arch,'PCWIN64') % WINDOWS
    ProjectionDataPath = 'C:\Users\rliu.MEDCTR\Desktop\L067\L067\full_DICOM-CT-PD';
    dictionaryFile =  'C:\Users\rliu.MEDCTR\Desktop\GeneralInfo\GeneralInfo\DICOM-CT-PD-dict_v8.txt';
elseif strcmp(arch, 'MACI64') % We haven't support MAC yet

end



Projection = readProjectionData(ProjectionDataPath, dictionaryFile);

xn = 128;
yn = 128;
zn = 128;
orgIdxX = (xn + 1.0) / 2.0; % MATLAB start with 1 therefore the center index should be 256.5
orgIdxY = (yn + 1.0) / 2.0;
orgIdxZ = 1;
ReconConf = genReconConf(Projection, xn, yn, zn, orgIdxX, orgIdxY, orgIdxZ);
useGPU = 0; % Use GPU

[ctimage, time] = analyticalRecon(Projection, ReconConf, useGPU);
ctimage1 = ctimage(:,:,:,1); % flying focal spot set 1
ctimage2 = ctimage(:,:,:,2); % flying focal spot set 2
img = ctimage1 + ctimage2;
imdisp(img);
