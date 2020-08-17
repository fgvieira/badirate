
# BadiRate

#### Estimating gene family turnover rates by likelihood-based methods
##### Copyright (C) 2012 Pablo Librado, Filipe Garrett Vieira, Julio Rozas and Universitat de Barcelona (UB)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).

## Synopsis 
This program estimates the gain, birth, death and innovation gene (or DNA element) family rates. It implements three stochastic population models: 

1. BDI: Birth-Death-Innovation model.

2. LI: Lambda-Innovation model, where Lambda is the birth-and-death parameter (equal birth and death rates are assumed)

3. GD: Gain-Death model, where gain parameter accounts for all gains regardless they are originated by DNA duplications or represent an innovation

4. BD: Birth and Death model

5. L: Lambda model, the birth-and-death parameter (equal birth and death rates are assumed). In addition, for each family or subfamily, it also infers the most likely ancestral content and the (sub)families unlikely evolving under the estimated stochastic process.

## Usage
The basic commands is:

    perl BadiRate.pl -tree NEWICK_FILE -fam FAMILY_FILE [options] > output.bd

or

    perl BadiRate.pl controlfile.bd > output.bd

where the options are specified in the control file (see manual for more details).

#### Available options:


    -anc		Reconstruct ancestral family sizes and minimum number of events
    -bmodel		Run local or free branch models
    -sizefile	Family Size File
    -unobs		Correct the likelihood for families absent in all extant species
    -h|help		Display this help
    -rmodel		Family turnover rates to be estimated
    -out		Set the output file
    -outlier	Report families escaping the estimated stochastic process
    -ep			Define the estimation procedure
    -print_ids	Display nodes ids in Newick format
    -priorfile	Prior File
    -n_max_int	Maximum allowed family size in internal phylogenetic nodes 
    -root_dist	Estimation method for the root a priori distribution
    -seed		Seed of the pseudo-random number generator
    -start_val	Starting values for the likelihood methods
    -treefile	Phylogenetic tree in Newick format
    -version	Report the BadiRate version

### Contact

For questions please contact the main developed at ![mail](https://bit.ly/2DYE86c)
