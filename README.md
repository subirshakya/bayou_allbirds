# bayou_allbirds

This uses the program bayou, a reverse jump MCMC method to implement an OU model on a tree for a given trait, for tarsus length data.

Tarsus data and tree used in analysis is part of the Avonet dataset that can be found here: https://figshare.com/ndownloader/files/34480865

Call R script in command line as:
```
Rscript ou_plot_cmd.R <num>
```
(here num refers to group number. I have divided the big bird dataset into 16 groups and included that as the first column of the dataset similar to 
https://onlinelibrary.wiley.com/doi/full/10.1002/evl3.267). 

The sixteen groups exactly as formatted in the csv file:
1. 1_Palaeognathae
2. 2_Galloanserae
3. 3_Phoenicopterimorpha
4. 4_Columbae         
5. 5_Strisores
6. 6_Gruiformes
7. 7_Aequorlitornithes
8. 8_Accipitriformes
9. 9_Coraciiformes
10. 10_Falconiformes
11. 11_Psittaciformes
12. 12_Suboscines
13. 13_Meliphagoidea
14. 14_Corvoidae
15. 15_Sylvoidea
16. 16_Muscicapoidea
17. 17_Nectaroidea
18. 18_Passeroidea   

Running rrun.sh <num> you can pass the number directly to the R script. Running multiple times shout create outfut files appended with iteration number automatically appened at the end. Or you can run a batch file and pass all values from 1 to 18 to run the whole thing with 18 parallel instances.
  
To get appropriate tuning parameters you can run Tuner.py which calls ou_test.R to get tuning parameters within the recommeneded 0.2 to 0.4 range. The starting group number and the starting values for alpha, beta, sigma, and theta can be updated in the python script.  
