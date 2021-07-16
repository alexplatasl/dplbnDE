dplbnDE
==========

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
run.DEbest
plot(run.DEbest)
```

To learn parameters of a custom structure, load a matrix of sizes edges x 2. Where 
columns represents direction (from-to) of edges. Like the followwing matrix:

``` r
#> my_structure
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

run.shade
plot(run.shade)
```

After the learning proccess. 
