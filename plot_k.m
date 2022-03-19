function plot_k2(kspace, segs, shots, fname)
    % kspace:   NCols, Nsegs, NPhases, Nshots, 3
    if nargin < 4
        fname = "untitled.png";
    fig = figure;
    for ii = 1:length(segs)
        scatter3( squeeze(kspace(:,segs(ii), 1, shots(ii), 1)), squeeze(kspace(:,segs(ii), 1, shots(ii), 2)), squeeze(kspace(:,segs(ii), 1, shots(ii), 3)) );
        axis([-pi pi -pi pi -pi pi]);
        p(ii) = "seg = "+num2str(segs(ii)) + "shot = "+num2str(shots(ii));
        legend(p);
        hold on;
        % drawnow;
    end
    saveas(fig, fname);
end