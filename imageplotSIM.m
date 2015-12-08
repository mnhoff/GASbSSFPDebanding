function [maxmg]=imageplotSIM(compdata,cyc,slrange,slind,numtis,th_pts,amat,bmat,noisSTD,flip,TR,T1,T2,GS,AS,GAS,CS,rowrep,CSgold,mmat,TRE_CS,TRE_GS,TRE_GAS,TRE_AS,var_GS,var_AS,var_GAS,varicompare,pr,cnt)
maxmg = 0.65*max(abs(compdata(:)));%for scaling
[nr,nc,nsl] = size(GS);
rr = 1;cc = 4;% for subplots
for sl = 1:nsl   
% for se = 1:4,figure(11);subplot(rr,cc,se);colormap(gray);imagesc(abs(compdata(:,:,se)), [0 maxmg]);
% title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);end 
% [~,h1]=suplabel(['Magnitude, Noise=',int2str(noisSTD)],'t');

%     for se = 1:4,figure(12);subplot(rr,cc,se);colormap(gray);imagesc(angle(compdata(:,:,se)),[-pi pi]);
%     title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);end
%     [~,h2]=suplabel(['Magnitude, Noise=',int2str(noisSTD)],'t');
%% just one
for se = 1:4
   figure;colormap(gray);imagesc(abs(compdata(:,:,se)),[0 maxmg]);
 axis image;set(gca,'xtick',[],'ytick',[]);
 title(['\Delta\theta=',int2str(cyc(se)*360),'\circ Original, {\alpha}=',...
    int2str(flip),'\circ, TR=',num2str(TR),'ms']);

    if numtis
        format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
            {[num2str(T1(1)),'/',num2str(T2(1))],...
             [num2str(T1(2)),'/',num2str(T2(2))],...
             [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
            [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
    else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
%         for shft = 1:10
%             text(2+16.1*(shft-1),137,'-\pi \rightarrow \pi','Color',[1 1 1]);
%         end
            set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
%         for shft = 1:5%-2pi to 2pi
%             text(4+32.2*(shft-1),137,'-2\pi\rightarrow2\pi','Color',[1 1 1],'FontSize',20,'FontWeight','bold');
%         end

    end
end

if pr
    COMP_str = ['Slice' num2str(slrange(sl)) '_' num2str(cyc(se)*2*180) 'PhaseCycle'];
    print(gcf,'-djpeg100',COMP_str)
            print(gcf,'-deps',COMP_str)
            saveas(gcf,COMP_str)
    close
end


%% plot M, for simdata
figure;colormap(gray);imagesc(mmat,[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
 title(['Base M, Noise=',int2str(noisSTD)],'FontSize',15);
if numtis
    format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
        {[num2str(T1(1)),'/',num2str(T2(1))],...
         [num2str(T1(2)),'/',num2str(T2(2))],...
         [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
        [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
set(gca,'XTick',1:8*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end

if pr
        M_str = 'M';
        print(gcf,'-djpeg100',M_str)
            print(gcf,'-deps',M_str)
            saveas(gcf,M_str)
        close
end

% %% plot CS gold, for simdata
figure;colormap(gray);imagesc(CSgold,[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
 title(['CS N = \inf, Noise=',int2str(noisSTD)],'FontSize',15);
if numtis
    format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
        {[num2str(T1(1)),'/',num2str(T2(1))],...
         [num2str(T1(2)),'/',num2str(T2(2))],...
         [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
        [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
else
    set(gca,'XTick',1:th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:th_pts:nc));
    set(gca,'YTick',nr/7:nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
    xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
    set(get(gca,'YLabel'),'Rotation',0.0)
end
if pr
    CSg_str = 'CS Gold Standard';
    print(gcf,'-djpeg100',CSg_str)
%             print(gcf,'-deps',CSg_str)
%             saveas(gcf,CSg_str)
    close
end


%% CS
figure;colormap(gray);imagesc(CS,[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
title(['CS, TRE = ',num2str(TRE_CS)],'FontSize',15);
if numtis
    format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
        {[num2str(T1(1)),'/',num2str(T2(1))],...
         [num2str(T1(2)),'/',num2str(T2(2))],...
         [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
        [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end

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
title(['GS, TRE = ',num2str(TRE_GS)],'FontSize',15);
if numtis
    format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
        {[num2str(T1(1)),'/',num2str(T2(1))],...
         [num2str(T1(2)),'/',num2str(T2(2))],...
         [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
        [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
set(gca,'XTick',1:8*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end
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
%     title(['GS Phase, Noise=',int2str(noisSTD)],'FontSize',15);
%     %         print -djpeg100 gs_phase.jpg
%% LGS = Linearized GS
% % % figure;colormap(gray);imagesc(abs(LGS(:,:,sl)),[0 maxmg]);
% % % axis image;set(gca,'xtick',[],'ytick',[]);
% % % title(['LGS, Noise=',int2str(noisSTD),' TRE = ',num2str(TRE_LGS)],'FontSize',15);
% % %     if numtis
% % %         format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
% % %             {[num2str(T1(1)),'/',num2str(T2(1))],...
% % %              [num2str(T1(2)),'/',num2str(T2(2))],...
% % %              [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
% % %             [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
% % % %         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
% % % %         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
% % %     else
% % % %         set(gca,'XTick',1:th_pts:nc)
% % % %         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
% % % %         set(gca,'YTick',nr/7:nr/7:nr)
% % % %         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
% % % %         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
% % % %         set(get(gca,'YLabel'),'Rotation',0.0)
% % %     set(gca,'XTick',1:8*th_pts:nc)
% % %     set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
% % %     set(gca,'YTick',nr/7:6*nr/7:nr)
% % %     set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
% % %     end
% % % if pr       
% % %      LGS_str = ['Slice' num2str(slrange(sl)) '_LGS'];
% % %      print(gcf,'-djpeg100',LGS_str)
% % % %          print(gcf,'-deps',LGS_str)
% % % %          saveas(gcf,LGS_str)
% % %      close
% % % end
%% Algebraics
%% no scaling - good if no singularities

figure;colormap(gray);imagesc(abs(AS(:,:,sl)),[0 maxmg]);
axis image;set(gca,'xtick',[],'ytick',[]);
title(['AS, TRE = ',num2str(TRE_AS)],'FontSize',15);
if numtis
    format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
        {[num2str(T1(1)),'/',num2str(T2(1))],...
         [num2str(T1(2)),'/',num2str(T2(2))],...
         [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
        [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
set(gca,'XTick',1:8*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end

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
figure;colormap(gray);imagesc(abs(GAS),[0 maxmg]);axis image;set(gca,'xtick',[],'ytick',[]);
title(['GAS, TRE = ',num2str(TRE_GAS)],'FontSize',15);
if numtis
    format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
        {[num2str(T1(1)),'/',num2str(T2(1))],...
         [num2str(T1(2)),'/',num2str(T2(2))],...
         [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
        [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)

set(gca,'XTick',1:8*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

end

if pr
    GAS_str = ['Slice' num2str(slrange(sl)) '_GAS'];
    print(gcf,'-djpeg100',GAS_str)
        print(gcf,'-deps',GAS_str)
        saveas(gcf,GAS_str)
    close
end


%% Variance maps
figure;colormap(gray);imagesc(var_GS,[0 max(var_AS(:))/2]);axis image;
title(['GS Variance, Sl', int2str(slind(1))]);
if ~numtis
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end
if pr
    GSvar_str = ['Slice' num2str(slrange(sl)) '_GS Variance'];
    print(gcf,'-djpeg100',GSvar_str)
        print(gcf,'-deps',GSvar_str)
        saveas(gcf,GSvar_str)
    close
end


figure;colormap(gray);imagesc(var_AS,[0 max(var_AS(:))/2]);axis image;
title(['AS Variance, Sl', int2str(slind(1))]);
if ~numtis
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end
if pr
    ASvar_str = ['Slice' num2str(slrange(sl)) '_AS Variance'];
    print(gcf,'-djpeg100',ASvar_str)
        print(gcf,'-deps',ASvar_str)
        saveas(gcf,ASvar_str)
    close
end

figure;colormap(gray);imagesc(var_GAS,[0 max(var_AS(:))/2]);axis image;
title(['GAS Variance, Sl', int2str(slind(1))]);
if ~numtis
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end
if pr
    GASvar_str = ['Slice' num2str(slrange(sl)) '_GAS Variance'];
    print(gcf,'-djpeg100',GASvar_str)
        print(gcf,'-deps',GASvar_str)
        saveas(gcf,GASvar_str)
    close
end

figure;colormap(gray);imagesc(varicompare);axis image;
title(['GAS Minimized Resistence in ',num2str(100*cnt/nr/nc),'% Pixels']);
if numtis
    format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
        {[num2str(T1(1)),'/',num2str(T2(1))],...
         [num2str(T1(2)),'/',num2str(T2(2))],...
         [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
        [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
    ylabel('T1/T2 (ms/ms)','Position',[-8 45],'FontSize',12);
    xlabel('\theta off-resonance','Position',[95 95],'FontSize',12);
else
    %this marks all the values
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
end
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
%     end
    
end