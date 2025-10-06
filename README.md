# Determine outgroup size
This script requires a .yaml file with a demography. The populations must be either named ingroup or out1 but there can be multiple outgroups (out1, out2, out3 and so on). 

The script simulated from a given demography. The user can specify how big the simulated chromosome is and how many chromosome to sample - in the examples below the chromosome size is 10 kb and we sample 1000 times - so a 10 Mb genome. There is uniform recombination rate.

Some examples:

```bash

python outgroup_efficientcy_github.py -outgroup_size=50 -chrom_size=10000 -iterations=1000 -demography=Demography_a.yaml -outfile=panel_a_simpledemography.pdf
python outgroup_efficientcy_github.py -outgroup_size=50 -chrom_size=10000 -iterations=1000 -demography=Demography_b.yaml -outfile=panel_b_simpledemography_3pop.pdf
python outgroup_efficientcy_github.py -outgroup_size=50 -chrom_size=10000 -iterations=1000 -demography=Demography_c.yaml -outfile=panel_c_exponentialgrowth.pdf
python outgroup_efficientcy_github.py -outgroup_size=50 -chrom_size=10000 -iterations=1000 -demography=Demography_d.yaml -outfile=panel_d_exponentialgrowth_3pop.pdf
python outgroup_efficientcy_github.py -outgroup_size=50 -chrom_size=10000 -iterations=1000 -demography=Demography_e.yaml -outfile=panel_e_ideal.pdf

```

The script will make a plot showing the demography you provided and how many of the non private SNPs were removed (closer to 100% the better)
