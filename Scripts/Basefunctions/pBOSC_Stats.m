 %% Evaluate model

 function stats = pBOSC_Stats(predictedepisodes, simulatedepisodes, powspctm, simdip, figflag)
%
 if ~exist('figflag')
    figflag = 0;
 end
%
diffepis = predictedepisodes + simulatedepisodes;
truepositive =  sum(sum(sum(diffepis==2))) ./ sum(sum(sum(simulatedepisodes))) * 100;
simfrx = find(sum(simulatedepisodes(simdip,:,:),3));
simpnts = find(simulatedepisodes(simdip,simfrx,:));


    % Figure: Model prediction
    if figflag
    cmp = colormap([1 1 1; [0 9 255]./255]);
    figure, imagesc(squeeze(predictedepisodes(simdip,:,:))), axis xy, axis off
    colormap(cmp)
    figure, imagesc(squeeze(simulatedepisodes(simdip,:,:))), axis xy, axis off
    colormap(cmp)
    end


% % Find the voxels with maximum number of the predicted episodes
[tmpval,leakvxidx] = sort(sum(predictedepisodes(:,simfrx,simpnts),3) / sum(simulatedepisodes(simdip,simfrx,simpnts)),'descend');
leakvxidx = leakvxidx(tmpval>.9); % Voxels with at least 90% of episodes found

% Now keep the voxel with maximum power
    [~,maxvox] = max(sum(powspctm(leakvxidx,simfrx,simpnts),3));
    maxvox = leakvxidx(maxvox);

% Stats simulation vs prediction
diffepis = predictedepisodes - simulatedepisodes;
falsenegative = sum(sum(sum(diffepis==-1))) ./ sum(sum(sum(simulatedepisodes))) * 100;
falsepositive = sum(sum(sum(diffepis==1)))  ./ sum(sum(sum(simulatedepisodes==0)))  * 100;
truenegative =  sum(sum(sum(predictedepisodes==0))) ./  sum(sum(sum(simulatedepisodes==0))) * 100;
precision = truepositive / (truepositive + falsepositive) * 100;
recall = truepositive / (truepositive + falsenegative) * 100;
f1score = 2 * (precision*recall) ./ (precision + recall);

% Alternative voxel
altvxepis = squeeze(predictedepisodes(maxvox,simfrx,:)) + squeeze(simulatedepisodes(simdip,simfrx,:));
truepositivealtvx =  sum(altvxepis==2) ./ sum( squeeze(simulatedepisodes(simdip,simfrx,:))) * 100;
altvxepis = squeeze(predictedepisodes(maxvox,simfrx,:)) - squeeze(simulatedepisodes(simdip,simfrx,:));
falsenegativealtvx =  sum(altvxepis==-1) ./ sum( squeeze(simulatedepisodes(simdip,simfrx,:))) * 100;
falsepositivevealtvx =  sum(altvxepis==1) ./  sum(sum(sum(simulatedepisodes==0))) * 100;
truenegativealtvx =  sum(sum(sum(predictedepisodes==0))) ./  sum(sum(sum(simulatedepisodes==0))) * 100;
precisionaltvx = truepositivealtvx / (truepositivealtvx + falsepositivevealtvx) * 100;
recallaltvx = truepositivealtvx / (truepositivealtvx + falsenegativealtvx) * 100;
f1scorealtvx = 2 * (precisionaltvx*recallaltvx) ./ (precisionaltvx + recallaltvx);

stats.simvoxel ={['True negative:']  [num2str(truenegative)];
                 ['False positive:'] [num2str(falsepositive)];
                 ['False negative:'] [num2str(falsenegative)];
                 ['True positive:']  [num2str(truepositive)];
                 ['Precision:']      [num2str(precision)];
                 ['Recall:']         [num2str(recall)]; 
                 ['F1 score:'] [num2str(f1score)]
                 ['Simulated voxel:'] [num2str(simdip)]};        


stats.altvoxel ={['True negative:']  [num2str(truenegativealtvx)];
                 ['False positive:'] [num2str(falsepositivevealtvx)];
                 ['False negative:'] [num2str(falsenegativealtvx)];
                 ['True positive:']  [num2str(truepositivealtvx)];
                 ['Precision:']      [num2str(precisionaltvx)];
                 ['Recall:']         [num2str(recallaltvx)]; 
                 ['F1 score:'] [num2str(f1scorealtvx)]
                 ['Alternative voxel:'] [num2str(maxvox)]};        
           
end


