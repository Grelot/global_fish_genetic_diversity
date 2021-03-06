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
## functions to calculate genetic diversity at species level 
## and cell level
## functions to calculate number of species, number of individuals per cell
##
#########################################################################
#Load the necessary library
using DataFrames, DataFramesMeta, StatsBase, CSV, Statistics

"""
Each line represents a cell (each cell is represented by a row)
with cell coordinates, mean, media and genetic diversity,
mean number of individuals, median of individuals, number of species,
mean number of individuals by species, median number of individuals by species
and standard deviation (sd)
"""
function calcGDmeanmednumberind(dats)
    #nb_indv
    nb_ind = by(unique(dats[[:cell, :species, :seq1]]),[:cell, :species]) do df
    tt = unique(df[:seq1])
    size(tt,1)+1
    end
    ind_mean = by(nb_ind, :cell, x -> mean(x[:x1]))
    ind_med = by(nb_ind, :cell, x -> median(x[:x1]))
    ind_std = by(nb_ind, :cell, x -> std(x[:x1]))
    #nb_species    
    richness = by(dats, :cell, totalspecies)
    #GD
    meanpaspec = by(dats, [:cell, :species], x -> mean(x[:num_per_bp]))
    GD_mean = by(meanpaspec, :cell, x -> mean(x[:x1]))
    GD_med = by(meanpaspec, :cell, x -> median(x[:x1]))
    #total
    total=DataFrame(site=GD_mean[:cell],
    GD_mean=GD_mean[:x1],
    GD_median=GD_med[:x1],
    nb_species=richness[:x1],
    nb_indv_mean=ind_mean[:x1],
    nb_indv_median=ind_med[:x1],
    nb_indv_sd=ind_std[:x1]
    )    
    total
end

"""
Calculate the GD (Genetic Diversity) value for each species for each cell

**Parameters**
* 'dats' : A master matrix data frame generated by the 'create_master_matrices()' function
"""
function calcGDspeciespersite(dats)
    nb_ind = by(unique(dats[[:cell, :species, :seq1]]),[:cell, :species]) do df
    tt = unique(df[:seq1])
    size(tt,1)+1
    end
    meanpaspec = by(dats, [:cell, :species], x -> mean(x[:num_per_bp]))
    total = DataFrame(site=nb_ind[:cell],
    species=nb_ind[:species],
    nb_indv=nb_ind[:x1],
    GDspecies=meanpaspec[:x1])
    total
end

"""
Calculate the GD (Genetic Diversity) value

**Parameters**
* 'dats' : A master matrix data frame generated by the 'create_master_matrices()' function
"""
function calcGD(dats)
    meanpaspec = by(dats, [:cell, :species], x -> mean(x[:num_per_bp]))
    GD = by(meanpaspec, :cell, x -> mean(x[:x1]))
    names!(GD, [:site, :GDval])
    GD, meanpaspec
end

"""
Calculate the weighted (by total # of sequence pairs) GD (Genetic Diversity) value

**Parameters**
* 'dats' : A master matrix data frame generated by the 'create_master_matrices()' function
"""
function calcGD_weighted(dats)
    meanpaspec = by(dats, [:cell, :species], x -> DataFrame(PIhat = mean(x[:num_per_bp]), numpairs = size(x,1)))
    GD = by(meanpaspec, :cell, x -> mean(x[:PIhat], WeightVec(x[:numpairs])))
    names!(GD, [:site, :GDval])
    GD, meanpaspec
end

"""
Calculate the total number of base pairs in a cell

**Parameters**
* 'dat' :A master matrix data frame generated by the 'create_master_matrices()' function
"""
function totalbasepairs(dat)
    seq1s = unique(dat[[:cell, :species, :seq1, :length_seq1]])
    seq2s = unique(dat[[:cell, :species, :seq2, :length_seq2]])
    seq=DataFrame(
    cell=vcat(seq1s[:cell],seq2s[:cell]),
    species=vcat(seq1s[:species],seq2s[:species]),
    seq=vcat(seq1s[:seq1],seq2s[:seq2]),
    length_seq=vcat(seq1s[:length_seq1],seq2s[:length_seq2])
    )
    seqUniq=unique(seq, [:species,:seq])
    sum(aggregate(seqUniq, [:cell,:species], sum)[:length_seq_sum])
end

"""
Calculate the total number of species in a cell

**Parameters**
* 'dat' :A master matrix data frame generated by the 'create_master_matrices()' function
"""
function totalspecies(dat)
    by(dat, :cell) do df
        sps = unique(df[:species])
        size(sps,1)
    end
end

"""
Calculate the richness and basepairs in a cell

**Parameters**
* 'dat' :A master matrix data frame generated by the 'create_master_matrices()' function
"""
function calc_cellvalues(dats)
    ret, full = calcGD(dats)
    ret[:richness] = by(dats, :cell, totalspecies)[:x1]
    ret[:basepairs] = by(dats, :cell, totalbasepairs)[:x1]
    ret, full
end


function cells_by_species(dats)
    by(dats,:species) do df
        unique(df[:cell])
    end
end



function ready_data(filename::AbstractString)
    ret = CSV.read(filename)
    ret = @where(ret, :seq1 .!= -1 )
    ret = @where(ret, :overlap .> 0.5)
    ret
end

