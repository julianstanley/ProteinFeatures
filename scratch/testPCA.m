% PCA plot of clusters.
blue = [0, 0, 255/255];
red = [255/255, 0, 0];
magenta = [255/255, 0, 255/255];
green = [0, 255/255, 0];
orange = [1, 0.4, 0];
colors = [blue; red; orange; magenta; green];

alnScoreMatrix = readtable("/home/julian/Dropbox (HMS)/Protein_Oxidation/FeatureEngineeringCode/Julian_Scripts/scratch/scores.csv");
alnScoreMatrix = table2array(alnScoreMatrix);
k = 5;
clusters = kmeans(alnScoreMatrix,k);

pcaCoeffs = pca(alnScoreMatrix,'NumComponents',2);
figure();
numClusters = max(max(clusters));
clusterSizes = [];
for k=1:numClusters
    clusterSizes(k,1) = size(find(clusters == k),1);
end
[Ks, Ki] = sort(clusterSizes,'descend');
for k=1:5
    inds = find(clusters == Ki(k));
    scatter(pcaCoeffs(inds,1),pcaCoeffs(inds,2),'filled','MarkerEdgeColor',colors(k,:),'MarkerFaceColor',colors(k,:));
    hold on;
    hull = convhull(pcaCoeffs(inds,1),pcaCoeffs(inds,2));
    kX = pcaCoeffs(inds,1);
    kY = pcaCoeffs(inds,2);
    plot(kX(hull),kY(hull),'-','Color',colors(k,:));
    hold on;
end
hold off;
