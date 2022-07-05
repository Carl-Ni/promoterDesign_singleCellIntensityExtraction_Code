function promoterData=GetPAO1promoterSequence(PAO1_DataStruct)
%   The function of this procedure is to obtain a sequence of bases upstream of each gene and not to determine whether it is a promoter.
%   Input parameter: PAO1_DataStruct is a structure containing the PAO1 genome sequence 'Fast' and annotation information 'Annotation'.
%   Output result: promoterData is also a structure, the structure contains:
%       Sequence：Sequence of bases upstream of the gene；
%       LocusTag：The locus tag for the corresponding gene；
%       fwdlength：Sequence contains the number of bases in the upstream gene；
%       Strand：Strand of the sequence；
%       fwdPrimerSequence：Sequence of forward primers；
%       fwdTm：Annealing temperature of forward  primers；
%       revPrimerSequence：Sequence of reserve primers；
%       revTm：Annealing temperature of reserve primers；
%       isgoodprimerpair：Is it a suitable primer pair；
%       AllLength：Length of promoter sequence;

%% Truncate the sequence fragment upstream of each gene without distinguishing whether it is a true promoter or not
promoterAllNums=length(PAO1_DataStruct.Annotation);
promoterDataNum=1;
for iPromoter=1:promoterAllNums
    if strcmp(PAO1_DataStruct.Annotation(iPromoter).Strand,'+')==true
        if iPromoter==1
            promoterData(promoterDataNum).LocusTag= PAO1_DataStruct.Annotation(iPromoter).LocusTag;
            promoterData(promoterDataNum).Sequence=PAO1_DataStruct.Fast(1:PAO1_DataStruct.Annotation(iPromoter).Location(1)+99);
            promoterData(promoterDataNum).Strand='+';
            promoterData(promoterDataNum).fwdlength=0;
            promoterDataNum=promoterDataNum+1;
        else
            if PAO1_DataStruct.Annotation(iPromoter).Location(1)-PAO1_DataStruct.Annotation(iPromoter-1).Location(2)>35
                promoterData(promoterDataNum).LocusTag= PAO1_DataStruct.Annotation(iPromoter).LocusTag;
                if PAO1_DataStruct.Annotation(iPromoter-1).Location(2)-PAO1_DataStruct.Annotation(iPromoter-1).Location(1)>1000
                    promoterData(promoterDataNum).Sequence=PAO1_DataStruct.Fast(PAO1_DataStruct.Annotation(iPromoter-1).Location(2)...
                        -999:PAO1_DataStruct.Annotation(iPromoter).Location(1)+99);
                    promoterData(promoterDataNum).fwdlength=1000;
                else
                    promoterData(promoterDataNum).Sequence=PAO1_DataStruct.Fast(PAO1_DataStruct.Annotation(iPromoter-1).Location(1)...
                        +1:PAO1_DataStruct.Annotation(iPromoter).Location(1)+99);
                    promoterData(promoterDataNum).fwdlength=PAO1_DataStruct.Annotation(iPromoter-1).Location(2)-PAO1_DataStruct.Annotation(iPromoter-1).Location(1);
                end
                promoterData(promoterDataNum).Strand='+';
                promoterDataNum=promoterDataNum+1;
            end
        end
    else
        if iPromoter==promoterAllNums
            promoterData(promoterDataNum).LocusTag= PAO1_DataStruct.Annotation(iPromoter).LocusTag;
            sequence=PAO1_DataStruct.Fast(PAO1_DataStruct.Annotation(iPromoter).Location(2)-99:end);
            promoterData(promoterDataNum).Sequence=seqrcomplement(sequence);
            promoterData(promoterDataNum).fwdlength=0;
            promoterData(promoterDataNum).Strand='-';
            promoterDataNum=promoterDataNum+1;
            clear sequence
        else
            if PAO1_DataStruct.Annotation(iPromoter+1).Location(1)-PAO1_DataStruct.Annotation(iPromoter).Location(2)>35
                promoterData(promoterDataNum).LocusTag= PAO1_DataStruct.Annotation(iPromoter).LocusTag;
                if PAO1_DataStruct.Annotation(iPromoter+1).Location(2)-PAO1_DataStruct.Annotation(iPromoter+1).Location(1)>1000
                    sequence=PAO1_DataStruct.Fast(PAO1_DataStruct.Annotation(iPromoter).Location(2)-99:PAO1_DataStruct.Annotation(iPromoter+1).Location(1)+999);
                    promoterData(promoterDataNum).fwdlength=1000;
                else
                    sequence=PAO1_DataStruct.Fast(PAO1_DataStruct.Annotation(iPromoter).Location(2)-99:PAO1_DataStruct.Annotation(iPromoter+1).Location(2)-1);
                    promoterData(promoterDataNum).fwdlength=PAO1_DataStruct.Annotation(iPromoter+1).Location(2)-PAO1_DataStruct.Annotation(iPromoter+1).Location(1);
                end
                promoterData(promoterDataNum).Sequence=seqrcomplement(sequence);
                promoterData(promoterDataNum).Strand='-';
                promoterDataNum=promoterDataNum+1;
                clear sequence
            end
        end
    end
end
%% Design of primers for the upstream sequences of the truncated genes
for kk=1:promoterDataNum-1
    [promoterData(kk).revPrimerSequence,promoterData(kk).revTm,promoterData(kk).isgoodprimerpair,...
        promoterData(kk).fwdPrimerSequence,promoterData(kk).fwdTm,promoterData(kk).AllLength]...
        =primerAutomated(promoterData(kk).Sequence,promoterData(kk).fwdlength);
    disp(kk)
end
end


function [revPrimerSequence,revTm,isgoodprimerpair,fwdPrimerSequence,fwdTm,Alllength]=primerAutomated(sequence,fwdlength)
%   Design of primers for the upstream sequences of the truncated genes
revPrimerSequence='a';
revTm=0;
isgoodprimerpair=0;
fwdPrimerSequence='a';
fwdTm=0;
Alllength=0;
compSequence=seqcomplement(sequence);
revSequence=s-eqreverse(compSequence);
%% All sequences of 18-27bp length were obtained
M=18:27;
for m=M
    fwdprimerSequence=sequence(1:fwdlength+m);
    revprimerSequence=revSequence(1:100+m);
    fwdindex=repmat((0:fwdlength)',1,m)+repmat(1:m,fwdlength+1,1);
    fwdprimerlist=fwdprimerSequence(fwdindex);
    revindex=repmat((0:100)',1,m)+repmat(1:m,100+1,1);
    revprimerlist=revprimerSequence(revindex);
    if m==M(1)
        for i=fwdlength+1:-1:1
            fwdprimerprops(i)=oligoprop(fwdprimerlist(i,:));
            fwdprimerEndNucleotide(i,1)=fwdprimerlist(i,end);
        end
        for ii=100+1:-1:1
            revprimerprops(ii)=oligoprop(revprimerlist(ii,:));
            revprimerEndNucleotide(ii,1)=revprimerlist(ii,end);
        end
    else
        for i=fwdlength+1:-1:1
            fwdprimerprops(end+1)=oligoprop(fwdprimerlist(i,:));
            fwdprimerEndNucleotide(end+1,1)=fwdprimerlist(i,end);
        end
        for ii=100+1:-1:1
            revprimerprops(end+1)=oligoprop(revprimerlist(ii,:));
            revprimerEndNucleotide(end+1,1)=revprimerlist(ii,end);
        end
    end
end
%% Screening for ineligible primers
%Filtering Primers Based on GC Content
fwdgc=[fwdprimerprops.GC]';
revgc=[revprimerprops.GC]';
bad_fwdprimers_gc=fwdgc<45|fwdgc>60;
bad_revprimers_gc=revgc<45|revgc>60;
%Filtering Primers Based on Their Melting Temperature
fwdtm=cell2mat({fwdprimerprops.Tm}');
revtm=cell2mat({revprimerprops.Tm}');
bad_fwdprimers_tm=fwdtm(:,5)< 55|fwdtm(:,5)>65;
bad_revprimers_tm=revtm(:,5)< 55|revtm(:,5)>65;
%Finding Primers Without a GC Clamp
bad_fwdprimers_clamp = lower(fwdprimerEndNucleotide(:,1)) == 'a' | lower(fwdprimerEndNucleotide(:,1)) == 't';
bad_revprimers_clamp = lower(revprimerEndNucleotide(:,1)) == 'a' | lower(revprimerEndNucleotide(:,1)) == 't';
%Finding Primers With Self-Dimerization and Hairpin Formation
bad_fwdprimers_dimers =~cellfun('isempty',{fwdprimerprops.Dimers}');
bad_fwdprimers_hairpin=~cellfun('isempty',{fwdprimerprops.Hairpins}');
bad_revprimers_dimers =~cellfun('isempty',{revprimerprops.Dimers}');
bad_revprimers_hairpin=~cellfun('isempty',{revprimerprops.Hairpins}');
bad_fwdprimers = [bad_fwdprimers_gc,bad_fwdprimers_tm,...
    bad_fwdprimers_dimers,bad_fwdprimers_hairpin,bad_fwdprimers_clamp];
bad_revprimers = [bad_revprimers_gc,bad_revprimers_tm,...
    bad_revprimers_dimers,bad_revprimers_hairpin,bad_revprimers_clamp];
good_fwdpos = find(all(~bad_fwdprimers,2));
good_fwdprop = fwdprimerprops(good_fwdpos);
good_revpos = find(all(~bad_revprimers,2));
good_revprop = revprimerprops(good_revpos);
if isempty(good_fwdpos)==false
    if isempty(good_revpos)==false
        for ifwd=1:length(good_fwdpos)
            fwdprimersequence=primersequence(good_fwdpos(ifwd),fwdlength+1,M(1),sequence);
            isNucleotideRepeats=NucleotideRepeats(fwdprimersequence);%Exclusion of sequences with 4 consecutive identical bases
            if isNucleotideRepeats==false
                for irev =1:length(good_revpos)
                    if abs(good_fwdprop(ifwd).Tm(5)-good_revprop(irev).Tm(5))<=5 %Look for primer pairs with a Tm difference of less than 5 degrees Celsius
                        revprimersequence=primersequence(good_revpos(irev),101,M(1),revSequence);
                        isNucleotideRepeats=NucleotideRepeats(revprimersequence);
                        if isNucleotideRepeats==false
                            revPrimerSequence=revprimersequence;
                            revTm=good_revprop(irev).Tm(5);
                            fwdPrimerSequence=fwdprimersequence;
                            fwdTm=good_fwdprop(ifwd).Tm(5);
                            isgoodprimerpair=1;
                            Alllength=length(sequence)-strfind(sequence,fwdPrimerSequence)-strfind(revSequence,revPrimerSequence)+2;%最终的PCR的序列
                            break
                        end
                    end
                end
                break
            end
        end
    end
end
end

function certainprimersequence=primersequence(pos,sequencelength,d0,sequence)
%Sequence of eligible primers obtained by positional information reduction
if mod(pos,sequencelength)~=0
    primerlength=d0+fix(pos/sequencelength);
    primeripos=mod(pos,sequencelength);
else
    primerlength=d0+fix(pos/sequencelength)-1;
    primeripos=sequencelength;
end
primerSequence=sequence(1:sequencelength-1+primerlength);
index=repmat((0:sequencelength-1)',1,primerlength)+repmat(1:primerlength,sequencelength,1);
primerlist=primerSequence(index);
certainprimersequence=primerlist(primeripos,:);
end

function isNucleotideRepeats=NucleotideRepeats(primer)
%Finding Primers With Nucleotide Repeats
repeats=regexpi(cellstr(primer),'a{4,}|c{4,}|g{4,}|t{4,}','ONCE');
isNucleotideRepeats=~cellfun('isempty',repeats);
end
