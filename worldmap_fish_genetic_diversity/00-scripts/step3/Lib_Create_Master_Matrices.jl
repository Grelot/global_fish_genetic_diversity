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
## Montpellier 2019
#########################################################################
##
## functions to create master data matrices 
## that are used to compute genetic diversity.
##
#########################################################################
# Load the necessary library and functions
include("Lib_Compare_Pairwise.jl")
using DataFrames
using DelimitedFiles
using CSV

"""
A function to create master data matrices that are used to compute genetic diversity, assess data quality and do sensitivity analyses.

**Parameters**
* 'foldername' : The name of of a folder containing the data files. There must be files of two types (file extension 'fasta' and file extension 'coords') with  the same filename, e.g. the species names (e.g. folder contents could be 'Bufo_bufo.fasta, Bufo_bufo.coords, Rana_arvalis.fasta, Rana_arvalis.coords', etc.). The .fasta files contain the sequences as an m x n integer matrix, where m is the number of sequences and n is the length of the alignments. The .coords files contain the geographic coordinates of the sequences, as an m x 2 floating point matrix with latitude in the first column and longitude in the second. 
"""
function create_master_matrices(foldername)
    cluster_name = basename(foldername)
    species_list = [x[1:(end-6)] for x in filter(st->occursin(".fasta",st), readdir(foldername))]      #identify unique file names ignoring the extension
    num_files = size(species_list,1)

    equalarea = DataFrame(species = String[], cell = String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])   # Pre-initialize the DataFrame to ensure correct element types
    latbands = DataFrame(species = String[], cell = String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])

    for iter_file = 1:num_files
        # read in the sequence data for one species
        file_name_seq = joinpath(foldername, species_list[iter_file] * ".fasta")        
        species_seqs = readdlm(file_name_seq, Int)
        length_seq = size(species_seqs,2)

        # read the coordinates for one species
        file_name_coords = joinpath(foldername, species_list[iter_file] * ".coords")
        species_coords = readdlm(file_name_coords)

        #EQUAL AREA GRID CELL
        file_name_equal = joinpath(foldername, species_list[iter_file] * ".equalareacoords")
        equ = [floor(Int, a) for a in readdlm(file_name_equal)[:,3]]
        equ = [String("$a") for a in equ]
        append!(equalarea, calcspecies(species_list[iter_file], species_seqs, equ))

        #LATITUDINAL BANDS
        ltbd = [string(10 * floor(Int,b/10)) for b in species_coords[:,1]]
        append!(latbands, calcspecies(species_list[iter_file], species_seqs, ltbd))
    end
    CSV.write("./07-master_matrices/$(cluster_name)_pairwise_equalarea.csv", equalarea)
    CSV.write("./07-master_matrices/$(cluster_name)_pairwise_latbands.csv", latbands)
end

"""
A function to calculate the summary statistic for all sites (e.g. grid cell or biome) for one species.

**Parameters**
* 'species':        A string with the name of the species
* 'species_seqs':   A matrix where the rows are aligned genetic sequences, and columns are loci. Basepairs must be coded as 1, 2, 3 or 4, or with a 0 signifying that the locus is absent from the alignment.
* 'sitenames':      A vector of strings with the names of the sites
"""
function calcspecies(species::String, species_seqs::Matrix{Int}, sitenames::Vector{String})
    res = DataFrame(species = String[], cell = String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])
    uniq_grid = unique(sitenames, dims=1)
    for iter_grid = 1:size(uniq_grid,1) # Go through each site in turn        
        seqs_grid = species_seqs[findall(x -> occursin(uniq_grid[iter_grid],x), sitenames),:]
        tot_muts = compare_pairwise(seqs_grid)
        lns = size(tot_muts, 1)
        tmp = DataFrame(species = fill(species, lns), cell = fill(uniq_grid[iter_grid, :][1], lns)) #Expand species and cell names to the length of the resulting DataFrame
        tot_muts_h = hcat(tmp, tot_muts)        
        append!(res, tot_muts_h)
    end

    #for sym in [:seq1, :seq2, :length_seq1, :length_seq2, :overlap, :commons, :num_per_bp]
    #   res[sym]=replace(res[sym], -1 => missing) #rep empty values with an NA identifier
    #end       
    res
end
