function recon2(ind, useTheta, discard)
    if nargin<2
        useTheta = 1;
    end
    if nargin<3
        discard = 1;
    end
include_path();
%% set the path and fname
fpath = "raw_data_15-3-22";
outpath = "tmp_recon_15-3-22";
% fpath = "raw_data_17-2-22";
% outpath = "recon_17-2-22";
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
% deadpts = deadptss(ind)
% issong = songorder(ind)
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

%% if retrospectively choose first nseg_retro
% if nargin<5
%     nseg_retro = Nsegs
% end
% k_traj = reshape(k_traj, NCols, Nsegs, Nshots, NPhases, 3);
% image = reshape(image, NCols, NCoils, Nsegs, Nshots, Navgs, NPhases);
% 
% Nsegs = 35;
% k_traj = reshape(k_traj(:, 2:end, :, :,:), NCols*Nsegs*Nshots, NPhases, 3);
% image = reshape(image(:,:,2:end, :,:,:), NCols, NCoils, Nsegs*Nshots, Navgs, NPhases);
% 
% NPhases = 1;
image = reshape(image, NCols, NCoils, Nsegs, Nshots, Navgs, NPhases);
plot_sig(image,8);
saveas(gcf,"untitled.png");
image = reshape(image(:,:,:, :,:,:), NCols, NCoils, Nsegs*Nshots, Navgs, NPhases);



%----------------------------------------------------
E = xfm_NUFFT([SpatialMtxSz, SpatialMtxSz, SpatialMtxSz, NPhases], [], [], k_traj);

%% reconstruct
for ii = 1:NCoils
    disp("reconstruct coil"+num2str(ii))
    tic
    % NCols, NCoils, NLines, Navgs, NPhases
    for jj = 1:NPhases
        kdata1(:,jj) = reshape(image(:,ii,:,1,jj), [], 1);
        kdata2(:,jj) = reshape(image(:,ii,:,2,jj), [], 1);
        kdata_diff(:,jj) = reshape(image(:,ii,:,2,jj) - image(:,ii,:,1,jj), [], 1);

        % grid1(:,:,:,jj,ii) = fftshift( reshape(E.st(jj).p' * squeeze(E.w(:,jj) .* kdata1(:,jj)), [mat,mat,mat].*2) );
        % grid2(:,:,:,jj,ii) = fftshift( reshape(E.st(jj).p' * squeeze(E.w(:,jj) .* kdata2(:,jj)), [mat,mat,mat].*2) );
    end
    recon1(:,:,:,:,ii) = reshape(E' *(E.w .* kdata1), mat, mat, mat, NPhases);
    recon2(:,:,:,:,ii) = reshape(E' *(E.w .* kdata2), mat, mat, mat, NPhases);
    recon_diff(:,:,:,:,ii) = reshape(E' *(E.w .* kdata_diff), mat, mat, mat, NPhases);

    % recon1(:,:,:,:,ii) = reshape(E.iter(kdata1, @pcg, 1e-4, 100, [1,1,1,0]), mat, mat, mat, NPhases);
    % recon2(:,:,:,:,ii) = reshape(E.iter(kdata2, @pcg, 1e-4, 100, [1,1,1,0]), mat, mat, mat, NPhases);
    % recon_diff(:,:,:,:,ii) = reshape(E.iter(kdata_diff, @pcg, 1e-4, 100, [1,1,1,0]), mat, mat, mat, NPhases);

    grid_diff(:,:,:,ii) = fftshift(reshape(E.st(1).p' * squeeze(E.w(:,1) .* kdata_diff(:,1)), [mat,mat,mat].*2));

    toc
end

% grid diff, weight
grid_weight = fftshift(reshape(E.st(1).p' * squeeze(E.w(:,1).^2), [mat,mat,mat].*2));

recon1 = squeeze(sum(abs(recon1).^2, 5).^0.5);
recon2 = squeeze(sum(abs(recon2).^2, 5).^0.5);
recon_diff = squeeze(sum(abs(recon_diff).^2, 5).^0.5);

grid_diff = squeeze(sum(abs(grid_diff).^2, 4).^0.5);

% grid1 = squeeze(sum(abs(grid1).^2, 5).^0.5);
% grid2 = squeeze(sum(abs(grid2).^2, 5).^0.5);
% grid1 = grid1./max(grid1(:));
% grid2 = grid2./max(grid2(:));
% size(grid1)
% size(grid2)

% ktmp = zeros(size(kdata1,1), 1);
% ktmp(1:NCols) = kdata1(1:NCols,1);
% tmp = fftshift(reshape(E.st(1).p' * squeeze(E.w(:,1) .* ktmp), [mat,mat,mat].*2));
% save_avw(abs(tmp), outpath + num2str(measID) + "_seg_" + num2str(nseg_retro) + "_tmp.nii.gz", 'd', [res/2, res/2, res/2]);

save_avw(abs(recon1), outpath + num2str(measID) + "_avg1.nii.gz", 'd', [res, res, res, 1]);
save_avw(abs(recon2), outpath + num2str(measID) + "_avg2.nii.gz", 'd', [res, res, res, 1]);

% save_avw(abs(grid1), outpath + num2str(measID) + "_grid1.nii.gz", 'd', [res/2, res/2, res/2, 1]);
% save_avw(abs(grid2), outpath + num2str(measID) + "_grid2.nii.gz", 'd', [res/2, res/2, res/2, 1]);

save_avw(abs(recon_diff), outpath + num2str(measID) + "_diff.nii.gz", 'd', [res, res, res, 1]);

save_avw(abs(grid_diff), outpath + num2str(measID) + "_grid_diff.nii.gz", 'd', [res/2, res/2, res/2]);
save_avw(abs(grid_weight), outpath + num2str(measID) + "_grid_weight.nii.gz", 'd', [res/2, res/2, res/2]);

disp("Finished")
end