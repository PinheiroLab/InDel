%% InDEL selection (v.0.1)
% By V. Pinheiro

% It uses the output of the inDEL assembly as input. As a minimum, it 
% requires a single cell array containing two columns, a peptide sequence
% as a string in column 1 and a frequency value in column 2.


for n = 1 : length(Freqs);
    Freqs{n} = horzcat('X', Freqs{n}, 'X');
end
Input = Freqs(:,1:2);
Output = Input;

% Creates input and output from selection matrices separate from the
% material generated in the assembly to avoid error, and it adds the 
% sequence context of the library

%% Simple Selection
% This is a very simple toy selection algorithm. It identifies sequences
% that contain the motif to be selected and multiplies their frequency
% by a chosen scalar

Sel_motif = {'XCR', 'XAF', 'YGDX', 'HGX'};
fold_enrich = [500, 100, 300, 100];

for a = 1: length(Sel_motif);
    for n = 1: length(Freqs);
        string_map = {};
        for b = 1: length(Freqs{n})-length(Sel_motif{a})+1;
            string_map {b,1} = Freqs{n}(b:b+length(Sel_motif{a})-1);
        end
        if ismember(Sel_motif{a}, string_map) == 1;
            Output{n,2} = Output{n,2} * fold_enrich(a) * rand();
        end
    end
end
% Routine to implement multiple concomitant slections

for n = 1: length(Freqs);
    string_map = {};
    for a = 1: length(Freqs{n})-length(Sel_motif)+1;
        string_map {a,1} = Freqs{n}(a:a+length(Sel_motif)-1);
    end
    if ismember(Sel_motif, string_map) == 1;
            Output{n,2} = Output{n,2} * fold_enrich * rand();
    end
end
% Breaks all available sequences in kmers the size of the selection motif
% being used and scans the whole library for its presence. Every time the
% motif is found, the frequency of the output is increased.

%% Subtracting Output from Input into selection

Enrichment = Output;


for n = 1: length(Input);
    [seq2, row] = ismember(Input(n), Output(:,1));
    if seq2 == 1;
        Enrichment {row, 2} = Enrichment {row, 2} - Input {n,2};
    else
        Input_intermediate = Input(n,:);
        Input_intermediate{1,2} = Input_intermediate{1,2} * -1; 
        Enrichment = cat(1, Enrichment, Input_intermediate);
    end
end
% This routine subtracts every sequence identified in Input that is still
% present in the Output. 

 
for n = length(Enrichment): -1: 1;
    if Enrichment{n, 2} == 0;
        Enrichment(n,:) = [];
    end
end
% This routine removes any sequence that is equally present in both Input
% and Output pools

varlist = {'a', 'b', 'Freqs', 'n', 'row', 'seq2', 'string_map', 'varlist'};
clear(varlist{:});
