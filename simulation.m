function simulation(ind, useTheta, discard)
    if nargin<2
        ind = 14
        useTheta = 1;
    end
    if nargin<3
        discard = 1;
    end
include_path();
%% set the path and fname
fpath = "raw_data_15-3-22/";
outpath = "sim_15-3-22/";

fpath = "/home/fs0/qijia/scratch/origin_data/cone_dev/phantom/" + fpath + "/"
outpath = "/home/fs0/qijia/scratch/origin_data/cone_dev/phantom/" + outpath + "/"
if ~isfolder(outpath)
    mkdir(outpath)
end

% gradnames = {'igurney_grad_64_20_9980', 'mjohnson_grad_64_200_9980', 'igurney_grad_64_20_9980', 'mjohnson_grad_64_200_9980', 'igurney_grad_128_20_10030', 'mjohnson_grad_128_200_10030', 'igurney_grad_128_20_5010','mjohnson_grad_128_200_5010'}
% deadpts = 1
codepath = pwd;
% measID = 220

%% load k-space data from Siemens .dat file
cd(fpath);
matchfile;  % load files matching gradient and measID, etc.

gradname = char(gradnames(ind))
measID = measIDs(ind)
if isstring(measID)
    measID = char(measID)
end
bwpixel = bwpixels(ind)

dead_ADC_pts = dead_ADC_ptss(ind)
twix_obj = mapVBVD(measID,'ignoreSeg',true,'removeOS',false);
cd(codepath);

if numel(twix_obj) > 1
    twix_obj = twix_obj{end};
end
dataSize = twix_obj.image.dataSize;

image=twix_obj.image(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
image = squeeze(image);
[NCols, NCoils, NLines, Navgs, NPhases] = size(image);
Nsegs = twix_obj.hdr.Config.NSeg;
Nshots = NLines / Nsegs;
KMtxSize = NCols;

%% load randnum
% randnum = load(fpath + "ktraj/randnum");
% randnum = randnum(1:Nsegs*NPhases*Nshots);
% Theta   = randnum / 32767.0 *2.0*pi;
if useTheta
randnum = twix_obj.image.iceParam(7,:);
randnum = reshape(randnum, Nsegs, NPhases, Navgs, Nshots);
randnum = reshape(randnum(:,:,1,:), [], 1);
Theta = randnum / 32767.0 *2.0*pi;
else
    Theta = zeros(Nsegs*NPhases*Nshots, 1);
end

%% set basic parameters
OS = twix_obj.hdr.Dicom.flReadoutOSFactor;
SpatialMtxSz = twix_obj.hdr.Config.BaseResolution;
mat = SpatialMtxSz;
fov = twix_obj.hdr.Config.RoFOV;
res = fov/SpatialMtxSz;
kmax = 5/res;

% bwpixel = 100;
Ts = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-6;
% bw_readout      = bwpixel * mat;                  % Hz
% Ts              = 1e3 / bw_readout / OS;               % ms, ADC sampling time
LEN             = mat * OS;
readout_time    = LEN * Ts;                           % ms
T = 10e-3;
grad_time       = ceil(readout_time / T) * T;
grad_pts = round(grad_time/T);
%% generate k-traj

if contains(gradname, "radial")
    base_k = zeros(grad_pts,3);
    base_k(:,3) = linspace(-kmax,kmax, grad_pts);
else
    base_g = load(fpath + "ktraj/" + gradname);
    % if length(base_g) < grad_pts
    %     deadpts = grad_pts - length(base_g);
        
        % disp("add dead points  "+num2str(deadpts));
    deadpts = ceil(dead_ADC_pts * Ts / 10e-3)

    base_g = [zeros(deadpts,3); base_g];
    % end
    base_k(:,1) = cumtrapz(squeeze(base_g(:,1))) .* 4.258 .* T;
    base_k(:,2) = cumtrapz(squeeze(base_g(:,2))) .* 4.258 .* T;
    base_k(:,3) = cumtrapz(squeeze(base_g(:,3))) .* 4.258 .* T;
end

base_k = [base_k(:,2), base_k(:,3), base_k(:,1)];   % Phase-Read-Slice coordinate system
% interpolation
base_k = calc_ADCpts(base_k, T, Ts, NCols);  % NCols x 3

kspace = zeros(NCols, Nsegs*NPhases*Nshots, 3);

% rotate base_k
GRCounter = twix_obj.image.iceParam(5,:);
GRCounter = reshape(GRCounter, Nsegs, NPhases, Navgs, Nshots);
GRCounter = reshape(GRCounter(:,:,1,:), [], 1);
[Azi, Polar] = GoldenMeans3D(GRCounter,true);
GrPRS = [sin(Azi).*sin(Polar), cos(Azi).*sin(Polar), cos(Polar)];
% [GrPRS, Azi, Polar, GRCounter] = gen_Grs(Nsegs, NPhases, Nshots, issong);
[GrPRS, GsPRS, GrRad, GsRad, R] = calc_slice(GrPRS, Theta);      % R [Nsegs*NPhases*Nshots, 3, 3]

for ii = 1: (Nsegs*NPhases*Nshots)
    kspace(:, ii, :) = (squeeze(R(ii,:,:)) * base_k')';
end

kspace = reshape(kspace, NCols, Nsegs, NPhases, Nshots, 3);
kspace = kspace./kmax.*pi;

%% discard dead ADC pts
if discard
    NCols = NCols - dead_ADC_pts;
    kspace = kspace(dead_ADC_pts+1:end, :, :, :, :);
    image = image(dead_ADC_pts+1:end, :, :, :, :);
    disp(dead_ADC_pts)
    size(kspace)
    size(image)
end

for ii = 1:NPhases
    k_traj(:,ii,:) = reshape(kspace(:,:,ii,:,:), [NCols*NLines, 1, 3]);
end

%----------------------------------------------------
tic
E = xfm_NUFFT([SpatialMtxSz, SpatialMtxSz, SpatialMtxSz, 1], [], [], k_traj(:,1,:));
toc


%% --------------------------------------------------------------------------------------------------------

%% calculate SNR, SNReff
under_factor = ceil(mat*mat*pi/2/NLines);
[w, kr, SNReff] = calc_dens(squeeze(k_traj(:,1,:)), mat, under_factor, pi);
densfig = plotMeanSDSlidingWindow(kr,abs(w),max(kr)/10,20,max(kr)*0.9, SNReff); 

[img, dims,scales,bpp,endian] = read_avw("/home/fs0/qijia/scratch/origin_data/cone_dev/phantom/sim_phantom_img/phantom_diff.nii.gz");
img = imresize3(img(:,:,:,1), [mat, mat, mat]);
kdata = E * img;
kdata = reshape(kdata, NCols, Nsegs, Nshots);
plot_kdata(kdata);
recon_img = reshape(E'* kdata(:), mat, mat, mat);
save_avw(abs(recon_img), outpath+num2str(measID)+"phantom_diff.nii.gz",'d', [res,res,res]);

% fname = 'radial_dn_2';
% load([dirname, fname, '.mat'], 'img', 'scales');
% img = imresize3(img, [mat, mat, mat]);
% recon_img = reshape(E.mtimes2(img), mat, mat, mat);
% save_avw(abs(recon_img), outpath+num2str(measID)+"_"+fname+".nii.gz",'d',[res,res,res]);

PSF = reshape(E' * ones(length(k_traj),1), [mat, mat, mat]);
[fwhm] = calc_fwhm(abs(PSF));
sidelobe_level = calc_sidelobe(abs(PSF));
save_avw(abs(PSF), outpath+num2str(measID)+"_PSF.nii.gz",'d',[res,res,res]);

% -----------------------------------
%%   compute SNR (pseudo replica)
noise_recons = add_noise(E, img);
[SNRreplica, SNRmap] = calc_SNR(1, noise_recons, img);
[SNRaliasing] = calc_SNR(3, noise_recons, img);
effSNRreplica = abs(SNRreplica/(prod(fwhm)));

saveas(densfig, outpath+num2str(measID)+"_DensityPlot.png");

save(outpath+num2str(measID)+"_metrics.mat", 'SNReff', 'fwhm', 'sidelobe_level', 'SNRreplica', 'effSNRreplica', 'SNRaliasing');

fprintf('SNReff = %.3f \n', abs(SNReff) );
fprintf('fwhm = %.2f, %.2f, %.2f\n', fwhm(1),fwhm(2),fwhm(3));
fprintf('sidelobe_level = %.4f, %.4f, %.4f\n', sidelobe_level(1),sidelobe_level(2),sidelobe_level(3));
fprintf('SNRreplica = %.2f\n', abs(SNRreplica));
fprintf('effSNRreplica = %.2f\n', abs(SNRreplica/(prod(fwhm))));
fprintf('SNRaliasing = %.2f\n', abs(SNRaliasing));

end