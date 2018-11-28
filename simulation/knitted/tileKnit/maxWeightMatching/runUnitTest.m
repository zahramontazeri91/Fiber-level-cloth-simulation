function passed = runUnitTest(edgeData,correct, inmaxcardinality)
if ~exist('inmaxcardinality')
    inmaxcardinality = false;
end

edgeData(:,1) = edgeData(:,1) +1; edgeData(:,2) = edgeData(:,2) +1; 
result = maxWeightMatching(edgeData, inmaxcardinality);
result = result-1;
result(find(result==-2)) = -1;
passed = all(result == correct);
% disp(result);
% disp(correct);
