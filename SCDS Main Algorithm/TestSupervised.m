% Train Supervised Standard Classifiers

%% Test Classifiers on Full Features set
Accuracy_KNN=[];
AccuracyStages_KNN=[];
Accuracy_RF=[];
AccuracyStages_RF=[];

for i=SubsetIds
    HypDat = Data.Hyp{i};
    FeatDat = FeatData{i};
    HypDat(HypDat==4)=3;
    HypDat(HypDat==5)=4;
    HypDat(HypDat==6)=5;
    HypDat(HypDat==0)=6;




    Labs_KNN = predict(Mod_KNN,FeatDat(HypDat<6,:));
    Labs_RF = predict(Mod_RF,FeatDat(HypDat<6,:));
    Labs_RF = str2double(Labs_RF);

    % plotconfusion(categorical(HypDat(HypDat<6)),categorical(Labs_KNN));
    cm = confusionmat(HypDat(HypDat<6),Labs_KNN);
    Accuracy_KNN = cat(1,Accuracy_KNN,sum(diag(cm))/sum(sum(cm))*100);
    AccuracyStages_KNN = cat(1,AccuracyStages_KNN,(diag(cm)./sum(cm,2)*100)');

    cm = confusionmat(HypDat(HypDat<6),Labs_RF);
    Accuracy_RF = cat(1,Accuracy_RF,sum(diag(cm))/sum(sum(cm))*100);
    AccuracyStages_RF = cat(1,AccuracyStages_RF,(diag(cm)./sum(cm,2)*100)');
end

Accuracy_KNN        = mean(Accuracy_KNN,1);
AccuracyStages_KNN  = mean(AccuracyStages_KNN,1);
Accuracy_RF         = mean(Accuracy_RF,1);
AccuracyStages_RF   = mean(AccuracyStages_RF,1);


%% MDF
FeatIDX=[];
for i=1:Cnum
    FeatIDX = [FeatIDX,FeatIDXMDF{i}];
end
FeatIDX = unique(FeatIDX);


Accuracy_KNN_MDF=[];
AccuracyStages_KNN_MDF=[];
Accuracy_RF_MDF=[];
AccuracyStages_RF_MDF=[];
for i=SubsetIds
    HypDat = Data.Hyp{i};
    FeatDat = FeatData{i}(:,FeatIDX);

    HypDat(HypDat==4)=3;
    HypDat(HypDat==5)=4;
    HypDat(HypDat==6)=5;
    HypDat(HypDat==0)=6;

    Labs_KNN_MDF = predict(Mod_KNN_MDF,FeatDat(HypDat<6,:));
    Labs_RF_MDF = predict(Mod_RF_MDF,FeatDat(HypDat<6,:));
    Labs_RF_MDF = str2double(Labs_RF_MDF);

    cm = confusionmat(HypDat(HypDat<6),Labs_KNN_MDF);
    Accuracy_KNN_MDF = cat(1,Accuracy_KNN_MDF,sum(diag(cm))/sum(sum(cm))*100);
    AccuracyStages_KNN_MDF = cat(1,AccuracyStages_KNN_MDF,(diag(cm)./sum(cm,2)*100)');

    cm = confusionmat(HypDat(HypDat<6),Labs_RF_MDF);
    Accuracy_RF_MDF = cat(1,Accuracy_RF_MDF,sum(diag(cm))/sum(sum(cm))*100);
    AccuracyStages_RF_MDF = cat(1,AccuracyStages_RF_MDF,(diag(cm)./sum(cm,2)*100)');
end

Accuracy_KNN_MDF        = mean(Accuracy_KNN_MDF,1);
AccuracyStages_KNN_MDF  = mean(AccuracyStages_KNN_MDF,1);
Accuracy_RF_MDF         = mean(Accuracy_RF_MDF,1);
AccuracyStages_RF_MDF   = mean(AccuracyStages_RF_MDF,1);