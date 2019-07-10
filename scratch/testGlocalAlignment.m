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

seq1 = "FVGGGVAFGP*PRSYDYTLPR";
seq2 = "GGGVAFGPKP*SYDYTLPRQV";

[scoreGlocal, alignmentGlocal] = nwalign(seq1,seq2,'Alphabet','AA','ScoringMatrix',hijack30,'Scale',1,'Glocal',true);
disp("Glocal = True (current parameter):")
disp(alignmentGlocal)
disp(scoreGlocal)

[scoreNGlocal, alignmentNGlocal] = nwalign(seq1,seq2,'Alphabet','AA','ScoringMatrix',hijack30,'Scale',1,'Glocal',false);
disp("Glocal = False (suggested parameter):")
disp(alignmentNGlocal)
disp(scoreNGlocal)