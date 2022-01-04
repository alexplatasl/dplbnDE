## `dplbnDE`: An R package for discriminative parameter learning of bayesian networks by Differential Evolution


Implements Differential Evolution (DE) to train parameters of Bayesian Networks 
(BN) for optimizing the Conditional Log-Likelihood (Discriminative Learning) 
instead of the log-likelihood (Generative Learning). Any given BN structure 
encodes assumptions about conditional independencies among the attributes and 
will result  in error if they do not hold in the data. Such an error includes 
the classification dimension. The main goal of Discriminative learning is 
minimize this type of error.

Install
=======

Make sure you have at least version 3.2.0 of R. You can get the current 
development version of `dplbnDE` from Github:

``` r
# install.packages('devtools')
devtools::install_github('alexplatasl/dplbnDE')
```

Example
=======

Load a data set and learn parameters of a bayesian network with custom structure 
or one learned by naive bayes, tree augmented naive Bayes using Chow-Liu’s algorithm 
or Hill-climbing.

``` r
library(dplbnDE)
data(car)
run.DEbest <- DEbest(NP=30, G=25, data = car, class.name = names(car)[7], crossover = "bin",
                mutation.pairs = 1, structure = "tan", F = 0.5, CR = 0.55,
                edgelist = NULL, verbose = 5)
#Gen:  5 	 CLL=  -1911.61 	 NP=  30 
#Gen:  10 	 CLL=  -1532.554 	 NP=  30 
#Gen:  15 	 CLL=  -1392.074 	 NP=  30 
#Gen:  20 	 CLL=  -1252.369 	 NP=  30 
#Gen:  25 	 CLL=  -1181.117 	 NP=  30 
                
run.DEbest
#Number of evaluations: 	 780 
#Final population size: 	 30 
#
#Summary results of fitness in final population: 
#
#Best CLL: 	 -1181.117 
#Worst CLL: 	 -1251.721 
#Median: 	 -1218.063 
#Std. Dev.: 	 17.96752 

#plot(run.DEbest)
```

To learn parameters of a custom structure, load a matrix of sizes edges x 2. Where 
columns represents direction (from-to) of edges. Like the following matrix:

``` r
my_structure
#     from       to        
#[1,] "class"    "buying"  
#[2,] "class"    "maint"   
#[3,] "class"    "doors"   
#[4,] "class"    "persons" 
#[5,] "class"    "lug_boot"
#[6,] "class"    "safety"  
#[7,] "maint"    "buying"  
#[8,] "lug_boot" "safety"

run.shade = lshade(NP=5, G=25, data = car, class.name = names(car)[7], c = 0.1,
             pB=0.05, edgelist = my_structure, verbose = 5)
#Gen:  5 	 CLL=  -1616.161 	 NP=  24 
#Gen:  10 	 CLL=  -1229.425 	 NP=  20 
#Gen:  15 	 CLL=  -1161.089 	 NP=  17 
#Gen:  20 	 CLL=  -1076.062 	 NP=  14 
#Gen:  25 	 CLL=  -1022.326 	 NP=  12 

run.shade
#Number of evaluations: 	 519 
#Final population size: 	 12 
#
#Summary results of fitness in final population: 
#
#Best CLL: 	 -1022.326 
#Worst CLL: 	 -1132.318 
#Median: 	 -1081.651 
#Std. Dev.: 	 39.99802 

#plot(run.shade)
```

After the learning process, returned bayesian networks can be analyzed with [`bnclassify`](https://cran.r-project.org/package=bnclassify) package.
