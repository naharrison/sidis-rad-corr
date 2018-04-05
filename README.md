# sidis-rad-corr
Radiative correction code for semi-inclusive DIS

### Requirements
* gfortran
* ...

### Instructions
* In haprad2-2014/Makefile, make sure CERNLIBS points to cernlib/x86_64_rhel7/2005/lib (see https://github.com/naharrison/cernlib-repo)
* To compile and run, do:
```
cd haprad2-2014
make
make test
./test.exe
```
* Output file is test.dat

### Comments
* The 7th argument of fhaprad is the minimum missing mass (e.g. ~1.2), not the average missing mass (e.g. ~1.596)
* Modified structure functions from Nick Markov
  * Pi+: H3_Updated = h3(default)* 12*(-0.03 + x); H4_Updated = 6*(2-x)*h4[default, bb = 5.1777]
  * Pi-: For H3 default values look reasonable.; H4_updated = 75*(2-x)*h4_default[bb = 2.077]
  * see: https://clasweb.jlab.org/rungroups/piprod/wiki/index.php/Haprad_calculations,_Nick,_07.16.2015#Comparison_of_structure_functions_.28as_a_function_of_Z.29
* Papers:
  * https://arxiv.org/abs/0711.4789
  * https://arxiv.org/abs/hep-ph/9903325
  
