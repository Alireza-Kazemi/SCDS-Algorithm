function DiscreteFlatPlot(X,Y,Labels)
%% invisible vertical line patch + dashed vertical lines
% prepare profile points, interleaving NaN between each pair
vnan = NaN(size(X)) ;
xp = reshape([X;vnan;X],1,[]); xp([1:2 end]) = [] ;
yp = reshape([Y;Y;vnan],1,[]); yp(end-2:end) = [] ;


% prepare the vertical lines, same method but we interleave the NaN at one
% element offset
xv = reshape([X;X;vnan],1,[]); xv([1:3 end]) = [] ;
yv = reshape([Y;vnan;Y],1,[]); yv([1:2 end-1:end]) = [] ;

% prepare colormap and color matrix (same method than above)
[~,~,colidx] = unique(Y) ;     
ncolor = length(Labels);           % Number of unique level
colormap(hsv(ncolor))           % assign a colormap with this number of color

% create the color matrix wich will be sent to the patch object
% same method of interleaving than for the X and Y coordinates
cd = reshape([colidx.';colidx.';vnan],1,[]); cd(end-2:end) = [] ;

% draw the patch (without vertical lines)
hp = patch(xp,yp,cd,'EdgeColor','flat','LineWidth',5) ;
% add the vertical dotted lines
hold on
hv = plot(xv,yv,':k') ;

% add a label centered colorbar
colorbar('Ticks',((1:ncolor)+.5)*ncolor/(ncolor+1),'TickLabels',Labels)
end