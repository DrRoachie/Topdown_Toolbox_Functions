function [clusters, clustersSUM, clustersMAX, clustersFreq]=ClusterMe(sv,fv)

%find clusters of threshold statistic test values ala E. Maris et al. (2007)
% sv is vector of thresholded stat values [either 0 or above 0]
% fv is vector of frequency values corresponding to sv
% output: clusters is a cell array with the thesholded stat-value clusters
% clustersSUM is the sum of each cell array
% clusterMAX is max sum
% clustersFreq is cell array w freq values corresponding to clusters

clear clusters
wrap = [0, sv, 0];
temp = diff(wrap ~= 0);
blockStart = find(temp == 1) + 1;
blockEnd = find(temp == -1);
clusters = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)),1:numel(blockStart), 'UniformOutput', false);


for i= 1:length(clusters)
clustersSUM(i) = sum(clusters{i});
end
clustersMAX = max(clustersSUM);

temp=find(sv==0);
fv(temp)=0;

clear clustersFreq
wrap = [0, fv, 0];
temp = diff(wrap ~= 0);
blockStart = find(temp == 1) + 1;
blockEnd = find(temp == -1);
clustersFreq = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)),1:numel(blockStart), 'UniformOutput', false);