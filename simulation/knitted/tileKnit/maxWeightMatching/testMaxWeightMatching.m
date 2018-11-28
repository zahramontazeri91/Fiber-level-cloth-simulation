% This code runs the unit tests on the maxWeightMatching function, almost
% all ported from Joris Van Rantwijk's unit tests.

result = maxWeightMatching([]);
result = maxWeightMatching([1 2 1]);

% Test with only one possible matching
edgeData = [2, 3, 10; 3, 4, 11];
result = maxWeightMatching(edgeData);

% Three edges
edgeData = [2 3 6; 3 4 12; 4 5 6];
result = maxWeightMatching(edgeData);

% Maximum cardinality
edgeData = [2 3 6; 3 4 12; 4 5 6];
result = maxWeightMatching(edgeData, false);

%Floating point weights
edgeData = [1 2 pi; 2 3 exp(1); 1 3 3.0; 1 4 sqrt(2)];
% edgeData(:,1) = edgeData(:,1) +1; edgeData(:,2) = edgeData(:,2) +1; 
runUnitTest(edgeData, [ -1, 4, 3, 2, 1 ])

% % Non-floating point version weights
edgeData = [1 2 4; 2 3 2; 1 3 3; 1 4 1];
runUnitTest(edgeData, [ -1, 4, 3, 2, 1 ], true)

% % Negative weights
edgeData = [1,2,2;1,3,-2;2,3,1;2,4,-1; 3,4,-6];
runUnitTest(edgeData, [ -1, 2, 1, -1, -1 ], false);
runUnitTest(edgeData, [ -1, 3, 4, 1, 2 ], true)


% % create S-blossom and use it for augmentation
edgeData = [1,2,8;1,3,9;2,3,10;3,4,7];
correct = [ -1, 2, 1, 4, 3 ];
runUnitTest(edgeData,correct)
edgeData = [1,2,8;1,3,9;2,3,10;3,4,7;1,6,5;4,5,6];
correct = [ -1, 6, 3, 2, 5, 4, 1 ];
runUnitTest(edgeData,correct)

% create S-blossom, relabel as T-blossom, use for augmentation
edgeData = [1,2,9;1,3,8;2,3,10;1,4,5;4,5,4;1,6,3];
correct = [-1, 6, 3, 2, 5, 4, 1 ];
runUnitTest(edgeData,correct)
edgeData = [1,2,9;1,3,8;2,3,10;1,4,5;4,5,3;1,6,4];
correct = [ -1, 6, 3, 2, 5, 4, 1 ];
runUnitTest(edgeData,correct)
edgeData = [1,2,9;1,3,8;2,3,10;1,4,5;4,5,3;3,6,4];
correct =  [ -1, 2, 1, 6, 5, 4, 3 ];
runUnitTest(edgeData,correct)

%create nested S-blossom, use for augmentation
edgeData = [1,2,9;1,3,9;2,3,10;2,4,8;3,5,8;4,5,10;5,6,6];
correct = [ -1, 3, 4, 1, 2, 6, 5 ];
runUnitTest(edgeData,correct)

%  create S-blossom, relabel as S, include in nested S-blossom
runUnitTest([1,2,10;1,7,10;2,3,12;3,4,20;3,5,20;4,5,25;5,6,10;6,7,10;7,8,8], [ -1, 2, 1, 4, 3, 6, 5, 8, 7 ])

%  create S-blossom, relabel as S, include in nested S-blossom
runUnitTest([1,2,10;1,7,10;2,3,12;3,4,20;3,5,20;4,5,25;5,6,10;6,7,10;7,8,8], [ -1, 2, 1, 4, 3, 6, 5, 8, 7 ])

% % create nested S-blossom, augment, expand recursively
runUnitTest([1,2,8;1,3,8;2,3,10;2,4,12;3,5,12;4,5,14;4,6,12;5,7,12;6,7,14;7,8,12],[ -1, 2, 1, 5, 6, 3, 4, 8, 7 ])

%create S-blossom, relabel as T, expand
runUnitTest([1,2,23;1,5,22;1,6,15;2,3,25;3,4,22;4,5,25;4,8,14;5,7,13],[ -1, 6, 3, 2, 8, 7, 1, 5, 4 ])

% %create nested S-blossom, relabel as T, expand
runUnitTest([1,2,19;1,3,20;1,8,8;2,3,25;2,4,18;3,5,18;4,5,13;4,7,7;5,6,7], [ -1, 8, 3, 2, 7, 6, 5, 4, 1 ])

% %create blossom, relabel as T in more than one way, expand, augment
runUnitTest([1,2,45;1,5,45;2,3,50;3,4,45;4,5,50;1,6,30;3,9,35;4,8,35;5,7,26;9,10,5], [ -1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9 ])


% again but slightly different
runUnitTest([1,2,45;1,5,45;2,3,50;3,4,45;4,5,50;1,6,30;3,9,35;4,8,26;5,7,40;9,10,5], [ -1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9 ])

% create blossom, relabel as T, expand such that a new least-slack S-to-free edge is produced, augment
runUnitTest([1,2,45;1,5,45;2,3,50;3,4,45;4,5,50;1,6,30;3,9,35;4,8,28;5,7,26;9,10,5], [ -1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9 ])

% create nested blossom, relabel as T in more than one way, expand outer blossom such that inner blossom ends up on an augmenting path
runUnitTest([1,2,45;1,7,45;2,3,50;3,4,45;4,5,95;4,6,94;5,6,94;6,7,50;1,8,30;3,11,35;5,9,36;7,10,26;11,12,5], [ -1, 8, 3, 2, 6, 9, 4, 10, 1, 5, 7, 12, 11 ])

% create nested S-blossom, relabel as S, expand recursively
runUnitTest([1,2,40;1,3,40;2,3,60;2,4,55;3,5,55;4,5,50;1,8,15;5,7,30;7,6,10;8,10,10;4,9,30], [ -1, 2, 1, 5, 9, 3, 7, 6, 10, 4, 8 ])

% Test case from Jan
%%
edgeData = [1 2 5; 1 3 6; 1 4 2; 2 3 4; 2 4 1;3 4 3];
result = maxWeightMatching(edgeData);
runUnitTest([1 2 5; 1 3 6; 1 4 2; 2 3 4; 2 4 1;3 4 3], [-1, 2, 1 ,4 ,3])

%%
edgeData = [1, 2, 47.2612; 1, 3, 46.9176; 2, 3, 49.3305; 1, 4, 44.7978; 2, 4, 49.1123; 2, 5, 51.1539; 4, 5, 50.5430; 2, 6, 48.2873; 3, 6, 47.7470; 4, 6, 46.8674; 5, 6, 48.8397]
result = maxWeightMatching(edgeData)
runUnitTest(edgeData, [-1 3     6     1     5     4     2])
