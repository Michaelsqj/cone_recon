function [Grs, Azi, Polar, GRCounter] = gen_Grs(Nsegs, NPhases, Nshots, issong)
    % generate readout gradient directions after rotation
    if issong
        GRCounter = reshape(0:(Nsegs*NPhases*Nshots-1), [], 1);
        GRCounter = reshape(GRCounter, Nshots, Nsegs, NPhases); 
        GRCounter = permute(GRCounter, [2,3,1]);    % Nsegs x Nphases x Nshots
    else
        GRCounter = reshape(0:(Nsegs*NPhases*Nshots-1), Nsegs, NPhases, Nshots);
    end
    
    GRCounter = GRCounter(:);
    [Azi, Polar] = GoldenMeans3D(GRCounter,true);

    Grs = [sin(Azi).*sin(Polar), cos(Azi).*sin(Polar), cos(Polar)];
end