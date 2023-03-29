# Omicron_BA.2

## Dataset

SARS-CoV-2 genomes and annotation information used in this study were downloaded from the GISAID database (https://www.gisaid.org) as of February 26, 2022. We collected 1,688,401 genomes which were 1) isolated from human samples and 2) annotated as Omicron (including BA.1, BA.1.1, BA.2 and BA.3 lineages). The number of undetermined nucleotides were counted for each genome, and 768,483 sequences were selected that 1) have certain sampling date, 2) were isolated from humans, 3) have less than 1000 undetermined nucleotides in its genome, and 4) have less than 10 undetermined nucleotides in the S protein region. For the BA.1 and BA.1.1 lineages, 286 genomes were used, which were sampled from 2021-08-12 to 2021-11-30. EPI_ISL_10023502 (BA.1) and EPI_ISL_10023526 (BA.1.1), which were both sampled from Republic of the Congo are the earliest samples in our data. For BA.2, 35 genomes were used, which were sampled from November 24, 2011 to December 10, 2021. The four genomes sampled from in November 2021 were detected in France (EPI_ISL_9796145, November 24, 2021), South Africa (EPI_ISL_8128463 and EPI_ISL_9679276, both November 27, 2021) and India (EPI_ISL_8693579, November 28, 2021). For BA.3, 28 genomes were used, which were sampled from November 24, 2011 to December 26, 2021. All of BA.3 variants sampled during November 2021 (in total 11 genomes) were isolated in South Africa. We also obtained non-omicron SARS-CoV-2 genomes by 1) two reference genomes [EPI_ISL_402125 (Wuhan-Hu-1, B lineage) and EPI_ISL_406862 (one of the earliest sequences carrying the S D614G mutation, B.1 lineage)]; 2) 20 randomly sampled genomes of each of the B.1.1.318 and B.1.1.519 lineages suggest by Majumdar et al.38 and Wang et al.39, respectively, and 3) randomly sampling five sequences for each month (from January 2021 to August 2021) for each 5 continents as Viana et al. conducted13. We excluded genomes that do not have PANGO categories or is assigned as recombinant of different PANGO lineages (i.e., lineage names starting from “X”). To further reduce the impact to recombination in data analysis, we had run recombination test using RDP4 software v4.101 (ref.40) multiple times. We excluded sequences that is involved in the recombination event, which has less than 3 sequences reported as recombinants. Finally, 349 Omicron and 275 non-Omicron genomes were used in this study. In total, 624 SARS-CoV-2 genomes used in this analysis were summarized in the following website: https://doi.org/10.55876/gis8.230110hk. 
