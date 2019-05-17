! gcc -std=c99 main.c; ./a.out

%%
clear
figure('Position',[400 400 1000 1000]);
axis tight square

X = load('output/xmat.txt');
Y = load('output/ymat.txt');
Z = load('output/zmat.txt');

ax = gca;
ax.NextPlot = 'replaceChildren';

surf(X,Y,Z)

% v = VideoWriter('surface_evolution_100Ma.avi');
% v.Quality = 95;
% open(v);
% % 
% % colormap gray;
% % for ii = 1:100
% %     Z = load(['output/zmat_',num2str(ii),'.txt']);
% %     pcolor(X,Y,gradient(Z)); shading flat; axis equal;
% %     pause(0.01);
% %     drawnow;
% %     writeVideo(v,getframe(gcf));
% % end
% % close(v);
% for ii = 1:3841
%     Z = load(['output/zmat_',num2str(ii),'.txt']);
%     ZZ(:,ii) = reshape(Z,1,[]);
% %     pcolor(X,Y,gradient(Z)); shading flat; axis equal;
% %     pause(0.01);
% %     drawnow;
% end
% disp('done');
% return;
% 
% subplot(1,2,1);
% % surf(X,Y,Z);
% 
% plot(Y(:,50), Z(:,50));
% subplot(1,2,2);
% imagesc(Z);
% axis equal;