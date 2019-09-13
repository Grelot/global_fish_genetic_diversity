##########################################################################
## Codes for the paper:
## Global determinants of freshwater and marine fish genetic diversity
## Authors :
## Stephanie Manel, Pierre-Edouard Guerin, David Mouillot,
## Simon Blanchet, Laure Velez, Camille Albouy, Loic Pellissier
## 
## Montpellier, 2017-2019
## Submited to Nature communications, 2019
##
##########################################################################
## This code is based on :
## An Anthropocene Map of Genetic Diversity
##
## Andreia Miraldo, Sen Li, Michael K. Borregaard,
## Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic,
## Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
##
## Submitted to Science, 2016
## Code in this file by Michael K. Borregaard
## modified by P.E. GUERIN
## Montpellier 2017
#########################################################################
##
## functions to compute the Genetic Diversity value from a set of sequences.
##
#########################################################################
# load the necessary library
using DataFrames
using DelimitedFiles

"""
A function to compute the GD (Genetic Diversity) value from a set of sequences. The algorithm compares all pairwise combinations of sequences and calculates the proportion of loci that differ between the pair. Returns a DataFrame with the identity of sequences, the lengths of each sequence, the overlap, and the computed pairwise divergence values.

**Parameters**
* 'species_seqs': A matrix where the rows are aligned genetic sequences, and columns are loci. Basepairs must be coded as 1, 2, 3 or 4, or with a 0 signifying that the locus is absent from the alignment.
"""
function compare_pairwise(species_seqs::Matrix{Int})
    num_seqs = size(species_seqs,1)  # count the number of sequences
    if(num_seqs < 2) # return an empty state if there is only one sequence
        return(DataFrame(seq1 = -1, seq2 = -1, length_seq1 = -1, length_seq2 = -1, overlap = -1., commons = -1, num_per_bp = -1.))
    end

    #Preinitialize containers
    seq1 = Int[]
    seq2 = Int[]
    length_seq1 = Int[]
    length_seq2 = Int[]
    overlap = Float64[]
    commons = Int[]
    num_per_bp = Float64[]
    for iter = 1:(num_seqs-1)
        for iter_sub = (iter + 1):num_seqs # For each pairwise combination of sequences
            pair_seqs = species_seqs[[iter, iter_sub], :]
            common_sites = minimum(pair_seqs,dims=1) .> 0 #Ensure that only sites existing in both sequences are compared (i.e. ignore loci missing from the alignment)
            if sum(common_sites) == 0.
                push!(num_per_bp, 0.)
            else
                num_segregating = sum(diff(pair_seqs[:, vec(common_sites)], dims=1) .!= 0)
                push!(num_per_bp, num_segregating/sum(common_sites))
            end
            push!(seq1, iter)
            push!(seq2, iter_sub)
            ls1 = sum(species_seqs[iter, :] .!= 0)
            ls2 = sum(species_seqs[iter_sub, :] .!= 0)
            push!(length_seq1, ls1)
            push!(length_seq2, ls2)
            push!(overlap, sum(common_sites)/max(ls1, ls2))  #Compute the overlap between the sequences
            push!(commons, sum(common_sites))
        end
    end
    DataFrame(seq1 = seq1, seq2 = seq2, length_seq1 = length_seq1, length_seq2 = length_seq2, overlap = overlap, commons = commons, num_per_bp = num_per_bp)
end
