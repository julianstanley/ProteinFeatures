function [scores, lengths] = pairwiseGlobalAlignment_1_13_2018_tempOnly(seqs1, seqs2)
% Manually paste in the sequences from tables or txt files, 1 seq per row
% after creating the following cell arrays...  Then this fxn can be called.
% seqs1 = cell(0);
% seqs2 = cell(0);
%

% Initialize return results.
%IDs = zeros(size(seqs2,1),1);
scores = zeros(size(seqs2,1),1);
lengths = zeros(size(seqs2,1),1);
%coverages = zeros(size(seqs2,1),1);
%matches = zeros(size(seqs2,1),1);

% Hijack the "*" of the blosum30 substitution matrix to force alignment
% with an "anchor point" or "fixpoint".
hijack30 = blosum30;
hijack30(24,24) = 1000;
% Modify "X" character in matrix to represent gap. Gap penalty = 8, and
% aligning an "X" with another "X" will be considered the same as aligning
% a gap with a gap, which normally wouldn't happen and will therefore be
% given a penalty of 0.
hijack30(23,:)=0;
hijack30(:,23)=0;
hijack30(23,23)=0;

%myCluster = parcluster('local');
%poolsize = myCluster.NumWorkers;
%parpool(poolsize);

% Loop through strands
parfor i=1:size(seqs2,1)

    % Get the two sequences to align.
    seq1 = char(seqs1{i,1});
    seq2 = char(seqs2{i,1});

    % Parse central residues:
    seq1center = seq1(ceil(length(seq1)/2):ceil(length(seq1)/2));
    seq2center = seq2(ceil(length(seq2)/2):ceil(length(seq2)/2));    

    % Replace central residues with "*" to force alignment of the central
    % residues.
    seq1 = strcat(seq1(1:ceil(length(seq1)/2)-1),'*',seq1(ceil(length(seq1)/2)+1:end));
    seq2 = strcat(seq2(1:ceil(length(seq2)/2)-1),'*',seq2(ceil(length(seq2)/2)+1:end));
    
	% Perform NW alignment
%	ID = eye(24);
%	[score, alignment] = nwalign(seq1,seq2,'Alphabet','AA','ScoringMatrix',ID,'Glocal',false,'GapOpen',realmax);
%	matches(i,1) = score;
	Length = min(length(seq2),length(seq1));
	lengths(i,1) = Length;
%	coverages(i,1) = size(alignment,2)/length(Length)*100;
%	IDs(i,1) = score/Length;
%	[score, alignment] = nwalign(seq1,seq2,'Alphabet','AA','ScoringMatrix','blosum30','Glocal',false,'GapOpen',realmax);
	[score, alignment] = nwalign(seq1,seq2,'Alphabet','AA','ScoringMatrix',hijack30,'Scale',0.2000,'Glocal',true);
    
    % Get score contribution from alignment of actual central residues.
	[scoreCenter, alignmentCenter] = nwalign(seq1center,seq2center,'Alphabet','AA','ScoringMatrix',hijack30,'Scale',0.2000,'Glocal',true);

    % Assemble the corrected score.
    scoreCorrected = score + scoreCenter - hijack30(24,24)*0.2;
	scores(i,1) = scoreCorrected;
end


