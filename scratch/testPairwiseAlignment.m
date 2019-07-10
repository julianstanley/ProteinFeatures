hijack30 = blosum30;
hijack30(24,24) = 1000;
hijack30(23,:)=0;
hijack30(:,23)=0;
hijack30(23,23)=0;


% Get the two sequences to align
    seq1 = char('AVDGFDIADAKTKNFISWAKQ');
    seq2 = char('AVDGFDIADAKTKNFISWAKQ');

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
    
%	coverages(i,1) = size(alignment,2)/length(Length)*100;
%	IDs(i,1) = score/Length;
%	[score, alignment] = nwalign(seq1,seq2,'Alphabet','AA','ScoringMatrix','blosum30','Glocal',false,'GapOpen',realmax);
	[score, alignment] = nwalign(seq1,seq2,'Alphabet','AA','ScoringMatrix',hijack30,'Scale',0.2000,'Glocal',true);
    
    % Get score contribution from alignment of actual central residues.
	[scoreCenter, alignmentCenter] = nwalign(seq1center,seq2center,'Alphabet','AA','ScoringMatrix',hijack30,'Scale',0.2000,'Glocal',false);

    % Assemble the corrected score.
    scoreCorrected = score + scoreCenter - hijack30(24,24)*0.2;
    
    disp(scoreCorrected)
    disp(alignment)
    disp(scoreCenter)
    

