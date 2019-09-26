# Hydrophobicity_Correlation

Hydrophobicity.py aims to calculate the similarity between a set of inputted/test sequences and a set of preestablished valid sequences. It does this through the following procedure:

1. Analyze the  preestablished valid sequences by areas of interest. Area of interest defined as consectutive hydropathy values that either all indicate a positive or negative hydrophobicity. With the default setting being 7, all strings of consecutive positive/negative hydropathy values with a length greater than 7 are stored as areas of interest.

2. Average out the stored areas of interest and create a range of hydropathy values for a relative position in the sequences. Assuming that all the sequences inputted are valid and are similar genetically, there should be a few regions of interest with a lot of different values for each.

3. Analyze the inputted/test sequences by areas of interest.

4. Determine whether or not the average for the area of interest for each test sequence lands within the range of hydropathy values for the same area of interest in the preestablished sequences. Using these determinations, the program creates a correlation value: (number of matched areas of interest with correct hydropathy values/total number of areas of interest)
Here the user can alter the results with +- hydropathy error value.

5. The program then filters out the sequences with the best correlation value and outputs them onto a FASTA file.

usage: Hydrophobicity.py [-h] [-w SIZE] [-e ERROR] [-n OUTPUT_SIZE]
                         [-v VIEW_MORE] [-o OVERLAP] [-g GRAPH]
                         input_reference input_test output

[-w Size] : modifies the size of interest for the windows (default is 7) -- potenitally needs to be modified in case sample size is too low

[-e Error] : modifies the error range when matching areas with the preestablished sequences (default is .1)

[-n Output_Size] : modifies the number of sequences desired from the output (default is 5)

[-v VIEW_MORE] : Lets you see the matched up areas

[-o Overlap] : Let's you modify the overlapping setting (default is 5 positions)

[-g GRAPH] : Graphs test sequences and a compilation of all reference sequences

input_reference: FASTA file format - a collection of sequences used to establish valid areas of interest and their hydropathy ranges to measure correlation.

input_test: FASTA file format - your test collection of sequences that are to be filtered out.

output: FASTA file format - the program writes these sequences onto the output.
