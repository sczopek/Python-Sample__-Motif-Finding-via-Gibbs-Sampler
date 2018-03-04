# DNA Motif Finding via Gibbs Sampler

This software demos the Gibbs Sampler algorithm by finding the Zinc Fingered GATA4 promoter motif in sample mouse DNA reads.  An overview for each file and the sample data is given, followed by some project notes including a getting started and installation guide.




### DemotGibbsSampler.py

This script demos the Gibbs Sampler Motif finding algorithm by finding the Zinc Fingered GATA4 promoter motif in sample mouse DNA reads using the using the methods in GibbsSampler.py.

### GibbsSampler.py

### Summary: 
This module implements the gibbs sampler algorithm, which is used to find common motifs in DNA sequences.  This probablistic search algorithm runs in polynomial time, which is an improvement over the brute force algorithm's exponential running time.
          
          
## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

**Python 2.7**

(Written in Python 2.7.12)

Ablility to run Python script from command line.

### Installing

Download and save all files to a local directory.

On the command line:
- Check for python -V 2.7
- Change to local directory containing the downloaded files.
- Run demo file to find Zinc GATA4 promoter in the mouse genome fragments.

```
python -V
CD /.../.../localDirWith Files
python DemoGibbsSampler.py
```


## Running the tests

This project was tested function by function as part of a CourseEra project using their web testing service.

The algorithm locates and prints each binding site.  You can compare the algorithm's answer to the real answer, by opening the solution file.  In that file the real motif binding site appears in capital letters.


## Versioning

v 0.1

## Authors

* **Scott Czopek** - *Initial work* - 1/22/17

## License

This project is free to copy and distribute.

## Acknowledgments

* I would like to thank my advisor Sergey Nuzdhin for giving me the chance to learn about genetics.
* Thanks to CourseEra I know more about motif finding.
