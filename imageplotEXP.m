function [maxmg]=imageplotEXP(compdata,cyc,slrange,slind,GS,AS,GAS,CS,var_GS,var_AS,var_GAS,varicompare,pr,cnt)
maxmg = 0.65*max(abs(compdata(:)));%for scaling, liver phantom
rr = 1;cc = 4;% for subplots
[nr,nc,nsl] = size(GS);
for sl = 1:nsl   
% for se = 1:4,figure(11);subplot(rr,cc,se);colormap(gray);imagesc(abs(compdata(:,:,se)), [0 maxmg]);
% title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);end 
% [~,h1]=suplabel(['Magnitude, Sl',int2str(sl+slind(1)-1)],'t');
% 
%     for se = 1:4,figure(12);subplot(rr,cc,se);colormap(gray);imagesc(angle(compdata(:,:,se)),[-pi pi]);
%     title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);end
%     [~,h2]=suplabel(['Phase, Sl',int2str(sl+slind(1)-1)],'t');
%% just one
for se = 1:4
   figure;colormap(gray);imagesc(abs(compdata(:,:,se,sl)),[0 maxmg]);
 axis image;set(gca,'xtick',[],'ytick',[]);
title(['\Delta\theta=',int2str(cyc(se)*360),'\circ Original, Sl',...
    int2str(sl+slind(1)-1)]);
end

if pr
    COMP_str = ['Slice' num2str(slrange(sl)) '_' num2str(cyc(se)*2*180) 'PhaseCycle'];
    print(gcf,'-djpeg100',COMP_str)
            print(gcf,'-deps',COMP_str)
            saveas(gcf,COMP_str)
    close
end

%% CS
figure;colormap(gray);imagesc(CS(:,:,sl),[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
title(['CS, Sl', int2str(sl+slind(1)-1)]);

if pr
        CS_str = ['Slice' num2str(slrange(sl)) '_CS'];
        print(gcf,'-djpeg100',CS_str)
            print(gcf,'-deps',CS_str)
            saveas(gcf,CS_str)
        close
end
%% GS
figure;colormap(gray);imagesc(abs(GS(:,:,sl)),[0 maxmg]);
axis image;set(gca,'xtick',[],'ytick',[]);
title(['GS Sl', int2str(sl+slind(1)-1)]);

if pr       
    GS_str = ['Slice' num2str(slrange(sl)) '_GS'];
    print(gcf,'-djpeg100',GS_str)
        print(gcf,'-deps',GS_str)
        saveas(gcf,GS_str)
         close
end
%     format_ticks(gca,{'-2\pi','-\pi','0','\pi','2\pi'},[],[-2*pi -pi 0 pi 2*pi]);

    %
%Phase
%     figure;colormap(gray);imagesc(angle(GS(:,:,sl)));axis image;set(gca,'xtick',[],'ytick',[]);
%     title(['GS Phase Sl', int2str(sl+slind(1)-1)],'FontSize',15);
%     %         print -djpeg100 gs_phase.jpg
% %% LGS = Linearized GS
% figure;colormap(gray);imagesc(abs(LGS(:,:,sl)),[0 maxmg]);
% axis image;set(gca,'xtick',[],'ytick',[]);
% title(['LGS Sl', int2str(sl+slind(1)-1)]);
% 
% if pr       
%      LGS_str = ['Slice' num2str(slrange(sl)) '_LGS'];
%      print(gcf,'-djpeg100',LGS_str)
% %          print(gcf,'-deps',LGS_str)
% %          saveas(gcf,LGS_str)
%      close
% end
%% Algebraics
%% no scaling - good if no singularities

figure;colormap(gray);imagesc(abs(AS(:,:,sl)),[0 maxmg]);
axis image;set(gca,'xtick',[],'ytick',[]);
title(['AS Sl', int2str(sl+slind(1)-1)]);

if pr
    AS_str = ['Slice' num2str(slrange(sl)) '_AS'];
    print(gcf,'-djpeg100',AS_str)
        print(gcf,'-deps',AS_str)
        saveas(gcf,AS_str)
    close
end
%
%Phase
%  figure;colormap(gray);imagesc(angle(AS(:,:,sl)));axis image;
%         title([int2str(npc),'cyc AS Phase Sl',int2str(sl+slind(1)-1)]);
%         print -djpeg100 as_phase.jpg
%% GAS
figure;colormap(gray);imagesc(abs(GAS(:,:,sl)),[0 maxmg]);axis image;set(gca,'xtick',[],'ytick',[]);
title(['GAS Sl', int2str(sl+slind(1)-1)]);

if pr
    GAS_str = ['Slice' num2str(slrange(sl)) '_GAS'];
    print(gcf,'-djpeg100',GAS_str)
        print(gcf,'-deps',GAS_str)
        saveas(gcf,GAS_str)
    close
end

%% Variance maps
figure;colormap(gray);imagesc(var_GS(:,:,sl),[0 max(var_AS(:))/2]);axis image;
title(['GS Variance, Sl', int2str(slind(1))]);

if pr
    GSvar_str = ['Slice' num2str(slrange(sl)) '_GS Variance'];
    print(gcf,'-djpeg100',GSvar_str)
        print(gcf,'-deps',GSvar_str)
        saveas(gcf,GSvar_str)
    close
end


figure;colormap(gray);imagesc(var_AS(:,:,sl),[0 max(var_AS(:))/2]);axis image;
title(['AS Variance, Sl', int2str(slind(1))]);

if pr
    ASvar_str = ['Slice' num2str(slrange(sl)) '_AS Variance'];
    print(gcf,'-djpeg100',ASvar_str)
        print(gcf,'-deps',ASvar_str)
        saveas(gcf,ASvar_str)
    close
end

figure;colormap(gray);imagesc(var_GAS(:,:,sl),[0 max(var_AS(:))/2]);axis image;
title(['GAS Variance, Sl', int2str(slind(1))]);

if pr
    GASvar_str = ['Slice' num2str(slrange(sl)) '_GAS Variance'];
    print(gcf,'-djpeg100',GASvar_str)
        print(gcf,'-deps',GASvar_str)
        saveas(gcf,GASvar_str)
    close
end

figure;colormap(gray);imagesc(varicompare(:,:,sl));axis image;
title(['GAS Min Variance in ',num2str(100*cnt(sl)/nr/nc),'% Pixels']);
if pr
    vartest_str = ['Slice' num2str(slrange(sl)) '_Variance Test'];
    print(gcf,'-djpeg100',vartest_str)
        print(gcf,'-deps',vartest_str)
        saveas(gcf,vartest_str)
    close
end
% 
%         %plot the number of pixels used in the variance calculation
%         figure;colormap(gray);imagesc(varpix,[0 9]);axis image;axis off;
%         title(['# Variance Pixels Employed, Sl', int2str(slind(1))]);
%                 print -djpeg100 numpix.jpg
%                 crop('numpix.jpg')
end   
