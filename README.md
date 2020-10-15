# SharkKinshipNetwork

README
Author, maintainer and contact.
Johann Mourier,
UMR MARBEC (IRD, Ifremer, Univ. Montpellier, CNRS), Sète, France. 
Email: johann.mourier@gmail.com

Description:

This repository includes all the R scripts to reproduce the analyses and plot the figures of the manuscript:
Mourier, J; Planes, S. 2020. Kinship does not predict the structure of a shark social network 

___________________________________________________________________________________________________________________________

Files

01_Main_GAI_SRI_Relatedness.R
Summary: This file contains the R Script to run the analyses and plot the figures of the manuscript.
02_functions.R
Summary: This file contains the R Script to run the custom functions required in the file 01_Main_GAI_SRI_Relatedness.R.

___________________________________________________________________________________________________________________________

Data
The data is available at the Dryad repository, where you can find the following files:


DataNetwork.csv

This file contains all observations of individual sharks (rows).
1.	ID contains the identification of all individuals observed in the groups;
2.	Group is the identification number of the groups;
3.	Date is the date at which the observation was made;
4.	Year is the year at which the observation was made;
5.	Month is the month number at which the observation was made;
6.	Period contains information on the period of observation (M = mating period, N = non-mating period);
7.	Site is the site ID where the observation was made;
8.	WeekNB is the week of the year at which the observation was made;
9.	Spweek is the continuous week number of these observations;
10.	SPMonth is the continuous month number of these observations
11.	Genetic contains information on whether genetic sample was taken for the individual.


DataBTRS_15sightings.csv

This file contains five variables of 156 groups (rows).
1.	date is the date at which the group was observed;
2.	Year is the year of the observed group;
3.	Month is the month of the observed group;
4.	Period contains information on the period of observation (M = mating period, N = non-mating period);
5.	Site is the site ID where the observation was made;
6.	Lat is the latitude of the site
7.	Long is the latitude of the site
8.	SP is the continuous week number of these observations
9.	IDs contains the identification of all individuals observed in the groups.


Attribute.csv

This file contains five covariables of 106 photo-identified individuals.
1.	ID contains the identity of all 106 photo-identified individuals;
2.	Sex contains the sex of each individual, where M refer to male individuals and F to females;
3.	Size contains the size (in cm) for each individual;
4.	Maturity contains the maturity status for each individual (MATURE = sexually mature, IMMATURE = sexually immature);
5.	DNA contains information on whether genetic sample was taken for the individual (Y = sampled, N = unsampled).

No.sightings.csv

This file contains three covariables of 106 photo-identified individuals in the overall study (without threshold).
1.	ID contains the identity of all 106 photo-identified individuals;
2.	Sightings contains the number of times the individual was sighted;
3.	Context indicates in which context the number of sightings corresponds to (All = overall study, Mating = mating period, NonMating = non-mating period);


MatrixGeneticRelatedness_Trioml.csv

This file contains the pairwise genetic relatedness (TrioML, Wang 2007) among 83 individuals.


SpatialProfile.csv

This file contains for each individual and each site an encounter rate (i.e., no. sightings of individual at site, divided by no. sampling occasions at site). 


CountSightings.csv

This file contains for each individual and each site the number of sightings. 


TempData.csv

This file contains five covariables to calculate temporal overlap.
1.	ID contains the identity of all 106 photo-identified individuals;
2.	Group is the identification number of the groups;
3.	Date is the date at which the observation was made;
4.	Site is the site ID where the observation was made;
5.	SPMonth is the continuous month number of these observations


Genotype(16markers).txt

This file contains the genotypes of 83 individuals at 16 microsatellite markers.


