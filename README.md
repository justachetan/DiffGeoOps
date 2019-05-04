# DiffGeoOps

This repository contains a python implementation of the paper:

>[Discrete Differential-Geometry Operators for Triangulated 2-Manifolds. *Mark Meyer*, *Mathieu Desbrun*, *Peter Schröder* and *Alan H. Barr*. *VisMath 2002*](http://www.multires.caltech.edu/pubs/diffGeoOps.pdf)

This module implements the discrete versions of three differential geometry operators that have been discussed in the paper. They are:

- Mean Curvature
- Gaussian Curvature
- Principal Curvature

## Usage

```
$ python3 DiffGeoOps.py -h
usage: DiffGeoOps.py [-h] [--calc] [--op OP] [-o O] [--mesh MESH]
                     [--title TITLE] [--mode MODE]
                     i

First, use --calc to generate files for containing value of the
operator and then plot the operatore using --mesh.

positional arguments:
  i              Path to input file

optional arguments:
  -h, --help     show this help message and exit
  --calc         Flag for calculation mode
  --op OP        Sets the operation
                 1 - Mean Curvature
                 2 - Gaussian Curvature
                 3 - Principal Curvatures
  -o O           Path to output file
  --mesh MESH    Required for plotting the figure
  --title TITLE  Title of plot
  --mode MODE    Tells how to handle obtuse triangles:
                 0: Voronoi Mode
                 1: Heron mode
                 2: Modified Heron Mode
```


First you need to provide a 3D mesh as an input to the code. The mesh should be in Object File Format (.off).
Use mode 2 for now until further notice. Then, use the command without the `--calc` flag for plotting by providing apprpriate mesh and value files.

### Example Usage

```
$ # Calculating Gaussian Curvature of a Torus mesh
$ python3 DiffGeoOps.py --calc --op 2 --mode 2 torus.off
$ # Plotting the generated values. The file generated will be <MESH-NAME>_<OP>.npy
$ python3 DiffGeoOps.py --mesh torus.off torus_KG.npy
```

## Some results
- Plot of Gaussian curvature for Torus
![Gaussian curvature for Torus](img/torus_KG.png?raw=true "Gaussian curvature of Torus")

<br> <br>
- Plot of Mean curvature for the Mother-son mesh
![Mean curvature for Mother-son mesh](img/mother.png "Mean curvature plot for Mother-son mesh")


## License 

Copyright (c) 2019 Aditya Chetan

For license information, see [LICENSE](LICENSE) or http://mit-license.org


- - -

This code was written as a part of my independent study in Differential Geometry with [Dr. Kaushik Kalyanaraman](https://www.iiitd.ac.in/kaushik) at IIIT Delhi during Winter 2019 Semester. 

For bugs in the code, please write to: aditya16217 [at] iiitd [dot] ac [dot] in


