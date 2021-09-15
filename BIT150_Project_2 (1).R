setwd('~/Desktop')
allele_freq <- read.csv(file = 'pop_allele_freq_gnomAD.csv') #import allel freq csv into R
View(allele_freq)
allele_freq_filter <- allele_freq[ -c(2:8,12,14,15) ] #filter out unneeded columnds from allele freq
allele_freq_filter <- transform(allele_freq_filter, AC_Afr_Am = Allele.Count.African.African.American / Allele.Number.African.African.American)
allele_freq_filter <- transform(allele_freq_filter, AC_Lat_Am = Allele.Count.Latino.Admixed.American / Allele.Number.Latino.Admixed.American)
allele_freq_filter <- transform(allele_freq_filter, AC_Ashkenazi_Jewish = Allele.Count.Ashkenazi.Jewish / Allele.Number.Ashkenazi.Jewish)
allele_freq_filter <- transform(allele_freq_filter, AC_Eas_As = Allele.Count.East.Asian / Allele.Number.East.Asian)
allele_freq_filter <- transform(allele_freq_filter, AC_Eur_Fin = Allele.Count.European..Finnish. / Allele.Number.European..Finnish.)
allele_freq_filter <- transform(allele_freq_filter, AC_Eur_non_Fin = Allele.Count.European..non.Finnish. / Allele.Number.European..non.Finnish.)
allele_freq_filter <- transform(allele_freq_filter, AC_South_As = Allele.Count.South.Asian / Allele.Number.South.Asian)
# above are calculations for allele frequencies for each population
allele_freq_filter <- allele_freq_filter[ -c(6:37) ] #filter out unneeded columnds from allele freq
View(allele_freq_filter)

                                

