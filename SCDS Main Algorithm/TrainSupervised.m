%% Train Supervised Standard Classifiers

HypDat = [];
FeatDat = [];
for i=SubsetIds
    HypDat = cat(1,HypDat,Data.Hyp{i});
    FeatDat = cat(1,FeatDat,FeatData{i});
end
HypDat(HypDat==4)=3;
HypDat(HypDat==5)=4;
HypDat(HypDat==6)=5;
HypDat(HypDat==0)=6;


Mod_KNN = fitcknn(FeatDat(HypDat<6,:), HypDat(HypDat<6), ...
                    'Distance', 'Euclidean', ...
                    'Exponent', [], ...
                    'NumNeighbors', 10, ...
                    'DistanceWeight', 'SquaredInverse', ...
                    'Standardize', false);
Mod_RF = TreeBagger(30,FeatDat(HypDat<6,:), HypDat(HypDat<6),'OOBPrediction','On',...
    'Method','classification');

%% MDF
FeatIDX=[];
for i=1:Cnum
    FeatIDX = [FeatIDX,FeatIDXMDF{i}];
end
FeatIDX = unique(FeatIDX);

HypDat = [];
FeatDat = [];
for i=SubsetIds
    HypDat = cat(1,HypDat,Data.Hyp{i});
    FeatDat = cat(1,FeatDat,FeatData{i}(:,FeatIDX));
end
HypDat(HypDat==4)=3;
HypDat(HypDat==5)=4;
HypDat(HypDat==6)=5;
HypDat(HypDat==0)=6;


Mod_KNN_MDF = fitcknn(FeatDat(HypDat<6,:), HypDat(HypDat<6), ...
                    'Distance', 'Euclidean', ...
                    'Exponent', [], ...
                    'NumNeighbors', 10, ...
                    'DistanceWeight', 'SquaredInverse', ...
                    'Standardize', false);
                
t = templateTree('Reproducible',true);
Mod_RF_MDF = TreeBagger(50,FeatDat(HypDat<6,:), HypDat(HypDat<6),'OOBPrediction','On',...
    'Method','classification');


                
