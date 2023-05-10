# hydroscatt



## Description
hydroscatt is a collection of python scripts to simulate the scattering properties of hydrometeors.

## Installation
The scripts require [pytmatrix](https://github.com/jleinonen/pytmatrix/) to run. All scripts use python3

To run the two-layer T-matrix code you need the f2py3 program from numpy:

- Set yourself into the f_2l_tmatrix directory

- Run the command:

 f2py3 -c tmatrix_py.f -m tm2l
 
 where tm2l is the name of the python module containing the fortran subroutines.

## Usage
This is a collection of python scripts and auxiliary functions. The scripts are meant to be self-explanatory. If something is not clear about their use or there are unexpected results please open an issue.

## Support
If you spot any bug or have suggestions about functionality please open an issue.

## Contributing
**Contributions are very welcomed!!!**

If you would like to contribute please open an issue to discuss your contribution.

Make a pull request so that the changes can be reviewed.

Please, make sure that the scripts are functional and they work as intended before uploading them. Comment them so that they are self-explanatory. Make sure that they comply with the [style guide for Python](https://peps.python.org/pep-0008/) by using a linter such as [pylint](https://www.pylint.org/).


## Authors and acknowledgment
This repository was created by Jordi Figueras.

It is based on IDL functions developed while working by MeteoSwiss. Please cite this paper if using hydroscatt in a publication:

- J. Figueras i Ventura, M. Boscacci, M. Gabella, U. Germann: Hydrometeor Classification at MeteoSwiss, in 8th European Conference on Radar in Meteorology and Hydrology ERAD2014, Garmisch-Partenkirchen, Germany, 1-5 Sept. 2014

The code uses PyTMatrix so if you use hydroscatt in a published work cite PyTMatrix accordingly.

The Fortran 2-layer T-matrix computations code was developed at Colorado State University with the advice of Larry Carey, Walt Petersen, and Bringi. The reference is:

- V. Bringi and T. Seliga, "Scattering from axisymmetric dielectrics or perfect conductors imbedded in an axisymmetric dielectric," in IEEE Transactions on Antennas and Propagation, vol. 25, no. 4, pp. 575-580, July 1977, doi: 10.1109/TAP.1977.1141642.

## License
BSD 3-Clause License

## Disclaimer
The software is still in a development stage. Please let us know if you would like to test it.

Open Radar cannot be held responsible for errors in the code or problems that could arise from its use.
