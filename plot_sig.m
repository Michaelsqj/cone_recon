function plot_sig(image,ind)

    % image NCols, NCoils, Nsegs, Nshots, Navgs, NPhases
    Nshots = size(image,4);    
    if nargin<1
       ind = 1:Nshots; 
    end
%     image = sqrt(sum(abs(image).^2, 2));
    size(image)
    figure;
    subplot(3,2,1);
    a1=squeeze(image(:,1,:,:,1,1));
    a1=reshape(a1(1:1:end,:,:),[],Nshots);
%     a1=abs(a1);
%     a1=abs(a1)/max(abs(a1(:)));
    plot(1:length(a1),abs(a1(:,ind)));
    ylim([0, max(abs(a1(:)))]);
    
    subplot(3,2,2);
    plot(1:length(a1),angle(a1(:,ind)));
    
    
    subplot(3,2,3);
    a2=squeeze(image(:,1,:,:,2,1));
    a2=reshape(a2(1:1:end,:,:),[],Nshots);
%     a2=abs(a2);
%     a2=abs(a2)/max(abs(a2(:)));
    plot(1:length(a2),abs(a2(:,ind)));
    ylim([0, max(abs(a2(:)))]);
    
    subplot(3,2,4);
    plot(1:length(a1),angle(a2(:,ind)));
    
    subplot(3,2,5);
    a=squeeze(image(:,1,:,:,1,1)-image(:,1,:,:,2,1));
    a=reshape(a(1:1:end,:,:),[],Nshots);
%     a=abs(a)/max(abs(a(:)));
%     a=abs(a);
    plot(1:length(a),abs(a2(:,ind))-abs(a1(:,ind)));
%     ylim([0, max(abs(a(:)))]);
    
    subplot(3,2,6);
    plot(1:length(a1),angle(a1(:,ind))-angle(a2(:,ind)));
    
%     figure;
%     for i=1:Nshots
%        plot(1:length(a), a(:,i));
%        title("shot = "+num2str(i));
%        ylim([0,1]);
%        xlim([-100,length(a)]);
%        drawnow;
%        pause(0.3);
%     end
end