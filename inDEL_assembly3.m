%% InDEL assembly prediction (v.0.2)
% By V. Pinheiro

%% Creating the assembly setup

prompt = ['How many assembly cycles'];
dlg_title = 'Assembly cycles';
num_lines = 1;
cycles = inputdlg(prompt, dlg_title, num_lines);
cycles = str2num(cycles{1});
% Asks user prompt for assembly cycle number

composition = zeros(20,cycles);

for n = 1: cycles;
    prompt = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    dlg_title = 'What is your desired composition for this cycle?';
    num_lines = 1;
    def = {'0', '0','0','0','0','0','0','0','0','0','0', '0','0','0','0','0','0','0','0','0'};
    composition_answer = inputdlg(prompt, dlg_title, num_lines, def);
    for a = 1: 20;
        composition(aa2int(prompt{a}),n) = str2num(composition_answer{a});
    end
end
% Asks user to choose ratio of amino acids for assembly cycle


composition_freq = zeros(20,cycles);
for n = 1: size(composition, 2);
    total = sum(composition(:,n));
    composition_freq (:,n) = composition (:,n)/total;
end
% Adjusting amino acid frequencies

figure(1);
imagesc(composition_freq, [0,1]);
set(gca,'xtick',[],'ytick',[]); %colormap bone;

%% Generating a population

prob_ligation = 0.5;
n_samples = 10000;
% Defining system variables

lib = zeros(n_samples,cycles);
for n = 1: n_samples;
    for x = 1: cycles;
        a = rand(1);
        if a <= prob_ligation;
            b = rand(1);
            for c = 1: 20;
       
                if b > 0;
                    lib(n,x) = c;
                    b = b - composition_freq(c,x);
                end
            end
        end
    end
end


% This for loop simulates the assembly method, generating a matrix (lib)
% containing the result of individual simulated assembly experiments. 
% Matrix uses the numerical codes for the amino acids in the cycle they are
% incorporated.


%% Collating the resulting libraries

Output = {};

for n = 1: n_samples;
    Clone = [];
    for x = 1: cycles;
        if lib (n,x) ~= 0;
            Clone = horzcat(Clone, int2aa(lib(n,x)));
        end
    end
    Output = vertcat(Output, Clone);
end

% This for loop converts the sparse simulation matrix into a cell
% containing the resulting 'obtained' library, which is the starting point
% of the analysis and similar to what would be obtained from a deep
% sequencing run.

%% Counting mutants

Freqs = tabulate(Output); 
% Make a table with the frequency of each variant

%% Creating a frequency and composition table

Max_column = ((cycles+1)/2)*cycles;
End_column = zeros(cycles, 1);
Start_column = zeros(cycles, 1);
Output_count = zeros (20, Max_column);
Output_frequency = zeros (20, Max_column);

% This creates the basic matrix where all loop lengths will be stored.
% Columns will be sequential so n=1 will be column 1, n=2 will be columns 2
% and 3, and so on. That will be the simplest route to making a master
% table that can be broken into smaller sets.


for n = 1: cycles;
    End_column(n) = ((n+1)/2)*n;
    Start_column(n) = End_column(n) - n +1;
end
% This generates a vector identifying the limits of each sublibrary

for n = 1: length(Output);
    a=0;
    for a = 1: length(Output{n});
        Column = ((length(Output{n})+1)/2)*length(Output{n}) - length(Output{n}) + a;
        Output_count(aa2int(Output{n}(a)),Column) = Output_count(aa2int(Output{n}(a)),Column) + 1;
    end
end

% Reads all assembly outputs separating them by length and counting the
% occurences of each amino acid at a given position of the library.

total = zeros(cycles,1);

% Creates a vector containing the sum of sequences for each sublibrary

for n = 1 : cycles;
    total(n) = sum(sum(Output_count(:,Start_column(n):End_column(n))));
    Output_frequency(:,Start_column(n):End_column(n)) = Output_count(:,Start_column(n):End_column(n))/total(n);
end

% Creates a frequency matrix for each sublibrary

total_frequencies = zeros(cycles,1);
for a = 1: length(total);
    if total(a) ~= 0;
        total_frequencies(a) = total(a)/sum(total);
    else
        total_frequencies(a) = 0;
    end
end

% Converts the total count into a frequency table

%% Histogram

figure(2);
bar(total_frequencies);

%% Amino acid positional probability per sublibrary

figure(3);

Gap = 0.5;
Left = zeros(cycles, 1);
Width = zeros(cycles, 1);

for n = 1: cycles;
    Left (n) = (Start_column(n) + (n*Gap))/(Max_column + (cycles+2)*Gap);
    Width (n) = (End_column(n) - Start_column (n) +1)/(Max_column + (cycles+2)*Gap); 
end

% This for loop defines the placement for each sublibrary graph

for b = 1: cycles;  
    if total(b) ~= 0;
        Max_freq = max(max(Output_frequency(:,Start_column(b):End_column(b))));
    else
        Max_freq = 1;
    end
    
    % This if routine normalises the display scale in each sublibrary
    % ensuring that the highest frequency incorporation is used to
    % normalise the displayed signal in each sublibrary
    
    axes('Position', [Left(b), 0.05, Width(b), 0.75]);     
    imagesc(Output_frequency(:,Start_column(b):End_column(b)), [0,Max_freq]);
    set(gca,'xtick',[],'ytick',[]); %colormap bone; 
    
    axes('Position', [Left(b), 0.85, Width(b), 0.05]);       
    imagesc(total_frequencies(b),[0,max(total_frequencies)]);
    set(gca,'xtick',[],'ytick',[]); 
end

% This generates a compound figure showing the frequency of the
% sublibraries and the positional composition of each sublibrary.

%% Amino acid positional probability taking the whole population into consideration

figure(4);

% This requires a new frequency matrix with positional distribution
% calculated for the entire population rather than per sublibrary

Output_frequency_all = Output_count/sum(total);
Max_freq2 = max(max(Output_frequency_all));

for b = 1: cycles;  
    
    
    axes('Position', [Left(b), 0.05, Width(b), 0.9]);     
    imagesc(Output_frequency_all(:,Start_column(b):End_column(b)), [0,Max_freq2]);
    set(gca,'xtick',[],'ytick',[]); %colormap bone; 
    
end
    
% This generates a compound figure showing the frequency of the
% sublibraries and the positional composition of each sublibrary based on
% the total frequency

%% Cleaning up unnecessary variables

varlist = {'Output', 'lib', 'a', 'b', 'c', 'Clone', 'Column', 'composition_answer', 'cycles', 'def', 'dlg_title', 'End_column', 'Gap', 'Left', 'Max_column', 'Max_freq', 'Max_freq2', 'n', 'num_lines', 'prompt', 'total', 'total_frequencies', 'Width', 'x', 'varlist'};
clear(varlist{:});


