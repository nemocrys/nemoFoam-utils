# nemocrys-foam

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10996331.svg)](https://doi.org/10.5281/zenodo.10996331)

OpenFOAM applications, utilities and boundary conditions applied in the NEMOCRYS project.

## Content

### Utilities for coupling of Elmer and OpenFOAM simulations

*  `scr/elmerFoamCoupling/utilities/elmerToFoamBC`  

   Maps boundary values from Elmer to an OpenFOAM mesh, including rotation to
3D depending on options provided in elmerToFoamDict. See https://github.com/nemocrys/elmerToFoamBC for more information.



* `scr/elmerFoamCoupling/solvers`  

    Modified versions of `buoyantBoussinesqPimpleFoam` and `buoyantBoussinesqSimpleFoam` with an additional Lorentz force and Joule heating distribution. The Lorentz force and Joule heat fields are defined in a separate csv-file and can be on a different mesh. The file `createFieldsFromElmer.H` handles the reading and mapping of the fields onto the OpenFOAM mesh. 
    The input fields can be 2D or 3D. Available mapping methods are nearest neighbor (`nearest`) and triangular linear interpolation (`interp1`).
For 2D input there is the option to map the field on an axisymmetric 3D simulation. This is done by finding the nearest neighbors of the target mesh projected on the 2D input plane.
See `tutorials/case3D` for an example. The available options are specified in the `elmerToFoamDict`file.
The Lorentz force and Joule heat are added as momentum and heat sources in the `Ueqn.H` and `Teqn.H` respectively.

*  `scr/elmerFoamCoupling/utilities/azimuthalAverage`  

    See `tutorials/azimuthalAverage` for an example setup. Utility to calculate the azimuthal average of the velocity on an axisymetric 3D mesh. The resulting velocity is mapped on a 2D wedge mesh and can be imported into Elmer. Both meshes should be available in the same case. 
    Given a 3D vector field in Cartesian coordinates $U_{3D} = (U_x, U_y, U_z)$, its projection onto the cylindrical coordinate system in the $x-y$ plane $U_{3D}^{proj} = (U_{r}, U_{\theta}, U_{z})$, can be computed as follows:
   
    $$U_r = \frac{U_x p_x + U_z p_z}{p_r}\   $$
    $$U_{\theta} = U_y\ $$
    $$U_z = 0 $$
   
    Where:
    - $U_x, U_y, U_z$ are the components of the 3D velocity field in Cartesian coordinates.
    - $p_x, p_y, p_z$ are the components of the position vector $\mathbf{p}_{3D}$ in Cartesian coordinates.

### boundary-conditions
* `shearStressMarangoni`  
See `tutorials/case3D` for an example.
Boundary condition with a Marangoni shear stress at a free surface due to a temperature gradient:
$$\nabla \mathbf{U} \cdot \mathbf{n} = -\frac{1}{\mu}\frac{\partial\gamma}{\partial T} \nabla T \cdot \mathbf{\tau}$$ 

## Compilation

To compile the solvers and boundary condition execute the `Allwmake`script. NOTE: The solvers and boundary conditions are configured to be put into the `FOAM_APPBIN` and `FOAM_LIBBIN`, respectively. This is required for seamless integration into [opencgs](https://github.com/nemocrys/opencgs).

## Referencing
If you use this code in your research, please cite our article (available with open access):

> Iason Tsiapkinis, Dadzis, K. (2024). nemocrys/nemoFoam-utils: Release of OpenFOAM applications, utilities and boundary conditions applied in the NEMOCRYS project. (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.10996331

## License

The code is provided under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.html).

## Acknowledgements

[This project](https://nemocrys.github.io/) has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 851768).

<img src="https://raw.githubusercontent.com/nemocrys/nemoblock/master/EU-ERC.png">

## Contribution

Any help to improve this package is very welcome!
