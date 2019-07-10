seqsTarget = cell(0);
seqsRef = cell(0);
% PASTE IN SEQUENCES

postProbs = [];

% PASTE IN POSTERIOR PROBABILITIES FOR ALL SEQS IN SAME ORDER AS IN
% ALIGNMENT SCORE MATRIX.

scores = zeros(size(seqsTarget,1),size(seqsRef,1));
lengths = zeros(size(seqsTarget,1),size(seqsRef,1));

tStart = tic;
for i=1:size(seqsRef,1)
    seqs2 = cell(size(seqsTarget,1));
    seqs2(:,1) = seqsRef(i,1);
    [score, length] = pairwiseGlobalAlignment_for_CS_clustering_7_17_2018(seqsTarget, seqs2);
    scores(:,i) = score;
    lengths(:,i) = length;
end
tElapsed = toc(tStart)


alnScoreMatrix = [];
alnScoreMatrix = scores;
postProbDiffsDups =[];
% OR PASTE IN SCORE MATRIX AND POSTERIOR PROBABILITY MATRIX.

% Clustergram of alignment scores.
CGobj = clustergram(alnScoreMatrix);
get(CGobj)

% Make matrix of positive/negative classes with same row/col order as
% clustergram for direct comparison/overlay.
% White = both row/col are positive
% Black = both row/col are negative
% Gray = 1 of row/col is positive and 1 is negative
classDiffMatrix = zeros(size(CGobj.RowLabels,1),size(CGobj.RowLabels,1));
for i=1:size(CGobj.RowLabels,1)

    tStart = tic;
    
    % Get the posterior probability corresponding to the same original row
    % as in the clustergram.
    rowPostProb = postProbs(str2num(CGobj.RowLabels{i,1}),1);    
    % Get the posterior probability corresponding to the same original
    % column as in the clustergram.
    for j=i:size(CGobj.ColumnLabels,2)
        colPostProb = postProbs(str2num(CGobj.ColumnLabels{1,j}),1);
        % Determine classDiffMatrix value.
        if ((rowPostProb == 1) && (colPostProb == 1))
            classDiffMatrix(i,j) = 1;
            classDiffMatrix(j,i) = 1;
        elseif ((rowPostProb == 0) && (colPostProb == 0))
            classDiffMatrix(i,j) = 0;
            classDiffMatrix(j,i) = 0;
        else
            classDiffMatrix(i,j) = 0.5;
            classDiffMatrix(j,i) = 0.5;
        end
    end
    i
    tElapsed = toc(tStart)
end

map = [0.0 0.0 0.0
    0.5 0.5 0.5
    1.0 1.0 1.0];
hmap = figure();
heatmap(flipud(classDiffMatrix));
colormap(hmap,map);

% Scatter plot of x = seq aln score, y = posterior probability difference.
alnScoreVector = [];
posProbDiffVector = [];
for i=1:size(alnScoreMatrix,1)
	alnScoreVector = vertcat(alnScoreVector,alnScoreMatrix((i+1):size(alnScoreMatrix,2),i));
	posProbDiffVector = vertcat(posProbDiffVector,postProbDiffsDups((i+1):size(postProbDiffsDups,2),i));
end

figure();
scatter(alnScoreVector(:,1),posProbDiffVector(:,1),50,'filled','o', ...
'MarkerFaceAlpha',0.05,'MarkerFaceColor','b');
hold on;
hold off;
[corrRHO,corrPVAL] = corr(alnScoreVector(:,1),posProbDiffVector(:,1),'type','Pearson');







% Cluster sequences based on alignment scores.
% Compute distances based on alignment scores.
D = pdist(alnScoreMatrix);
% Compute agglomerative hierarchical cluster tree from distance matrix.
tree = linkage(alnScoreMatrix,'average','euclidean');
leafOrder = optimalleaforder(tree,D);
figDendogram = figure();
[H,T,outperm] = dendrogram(tree,0,'reorder',leafOrder);
hAxis = get(H(1),'parent');

c = cophenet(tree,D);







% Find optimal cluster assignments by looping through all possible number
% of clusters and evaluating the ARI.
bestNumClusters = 0;
maxARIscore = -1*realmax;
ARIscores = [];
for n=1:size(alnScoreMatrix,1)
    n
	tStart = tic;
    
    % Cluster sequences by trimming the cluster tree to get 'n' clusters.
    T = cluster(tree,'maxclust',n);
    clusters = T;
    ARIscore = adjrand(postProbs,clusters);

    % Get unique cluster IDs.
    clusterIDs = unique(T);
%    % Loop through cluster assignments.
%    for i=1:size(clusterIDs,1)
%        % Identify next cluster number.
%        clustNum = clusterIDs(i,1);
%        % Get all indices in T that are in the same cluster.
%        inds = find(T == clustNum);
%        clustNumPosInds = size(find(postProbs(inds,1) > 0.5),1);
%        clustNumNegInds = size(find(postProbs(inds,1) < 0.5),1);
%        clustHamScore = ((clustNumPosInds^2)-clustNumPosInds)/2+((clustNumNegInds^2)-clustNumNegInds)/2-(clustNumPosInds*clustNumNegInds);
%        clustHamScore = ((clustNumPosInds^2)-clustNumPosInds)/2-(clustNumPosInds*clustNumNegInds);
%        HammingScore = HammingScore + clustHamScore;
%    end
    % Update best clustering stats.
    if (ARIscore > maxARIscore)
        maxARIscore = ARIscore;
        bestNumClusters = n;
    end
    
    % Store score for this cluster assignment.
    ARIscores(n,1) = ARIscore;
    
	tElapsed = toc(tStart)
end






% K-means clustering:
% Find optimal cluster assignments by looping through all possible number
% of clusters and evaluating the ARI.
ARIscores = [];
bestKclusters = [];
bestARIscore = 0;
for n=1:size(alnScoreMatrix,1)clusters
    n
	tStart = tic;
    
    % Cluster sequences by k-means.
    clusters = kmeans(alnScoreMatrix,n);
    ARIscore = adjrand(postProbs,clusters);
    
    % Store score for this cluster assignment.
    ARIscores(n,1) = ARIscore;
    
	if (ARIscore > bestARIscore)
        bestARIscore = ARIscore;
        bestKclusters = clusters;
    end

	tElapsed = toc(tStart)
end




% K-means clusters are variable even with same parameter settings. Repeat a
% defined number of trials and keep the best result.
ARIscores = [];
bestKclusters = [];
bestARIscore = 0;
k = 5;
trials = 10;
for n=1:trials
    n
	tStart = tic;
    
    % Cluster sequences by k-means.
    clusters = kmeans(alnScoreMatrix,k);
    ARIscore = adjrand(postProbs,clusters);
    
    % Store score for this cluster assignment.
    ARIscores(n,1) = ARIscore;
    
    if (ARIscore > bestARIscore)
        bestARIscore = ARIscore;
        bestKclusters = clusters;
    end
        
	tElapsed = toc(tStart)
end






% Testing stuff...

%rowLabels = [];
%for i=1:size(CGobj.RowLabels,1)
%   rowLabels(i,1) = str2num(CGobj.RowLabels{i,1});
%end
%sortedScoreMatrix = alnScoreMatrix(rowLabels,:);
%sortedScoreMatrix = sortedScoreMatrix(:,rowLabels);

[As, Ai] = sort(clusters);
[Bs, Bi] = sort(alnScoreMatrix(:, 1));
ABi(Ai) = Bi;

[a_sorted, a_order] = sort(clusters);
sortedScoreMatrix = alnScoreMatrix(a_order,:);
sortedScoreMatrix = sortedScoreMatrix(:,a_order);
figure();
HMobj = HeatMap(sortedScoreMatrix,'DisplayRange',3);
get(HMobj)




% Make matrix of positive/negative classes with same row/col order as
% clustergram for direct comparison/overlay.
% White = both row/col are positive
% Black = both row/col are negative
% Gray = 1 of row/col is positive and 1 is negative
classDiffMatrix = zeros(size(a_order,1),size(a_order,1));
for i=1:size(a_order,1)
    % Get the posterior probability corresponding to the same original row
    % as in the clustergram.
    rowPostProb = postProbs(a_order(i,1),1);
    % Get the posterior probability corresponding to the same original
    % column as in the clustergram.
    for j=i:size(a_order,1)
        colPostProb = postProbs(a_order(j,1),1);
        % Determine classDiffMatrix value.
        if ((rowPostProb == 1) && (colPostProb == 1))
            classDiffMatrix(i,j) = 1;
            classDiffMatrix(j,i) = 1;
        elseif ((rowPostProb == 0) && (colPostProb == 0))
            classDiffMatrix(i,j) = 0;
            classDiffMatrix(j,i) = 0;
        else
            classDiffMatrix(i,j) = 0.5;
            classDiffMatrix(j,i) = 0.5;
        end
    end
end

map = [0.0 0.0 0.0
    0.5 0.5 0.5
    1.0 1.0 1.0];
hmap = figure();
heatmap(flipud(classDiffMatrix));
colormap(hmap,map);


% PCA plot of clusters.
blue = [0, 0, 255/255];
red = [255/255, 0, 0];
magenta = [255/255, 0, 255/255];
green = [0, 255/255, 0];
orange = [1, 0.4, 0];
colors = [blue; red; orange; magenta; green];

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
