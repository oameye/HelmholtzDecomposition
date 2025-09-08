# HelmholtzDecomposition Package

A Mathematica package for performing analytical Helmholtz decomposition of vector fields in n dimensions, following the approach in Richters and Glotz, "Analytical Helmholtz Decomposition of n-Dimensional Vector Fields" ([10.1016/j.jmaa.2023.127138](https://doi.org/10.1016/j.jmaa.2023.127138)).


## Quick Start

```mathematica
(* Load the package *)
Needs["HelmholtzDecomposition`"];

(* Define your coordinate variables *)
x = Subscript[x, 1]; y = Subscript[x, 2];

(* Your vector field *)
f = {x^2 + y, x*y};
coords = {x, y};

(* Complete decomposition with output *)
PrintHelmholtz[f, coords]

(* Direct access to components *)
{gradient, rotational} = HelmholtzDecomposition[f, coords];
```

## Main Functions

One defines variables:
- `f`: Vector field as a list (e.g., `{x^2, y^2}`)
- `xvec`: Coordinate variables as a list (e.g., `{x, y}`)

The package as the following main functions:
- `HelmholtzDecomposition[f, xvec]`: Returns `{gradient, rotational}` components as a list.

- `HelmholtzGradient[f, xvec]`: Returns only the gradient component.

- `HelmholtzRotational[f, xvec]`: Returns only the rotational component.

- `PrintHelmholtz[f, xvec]`: Displays the complete Helmholtz decomposition as in the original code of Richters and Glotz.

## Citation and License

Any intellectual property rights of the original authors are fully respected. Please cite the original paper when using this code.
- [Zenodo repository](https://zenodo.org/records/7680297) 
- Erhard Glötzl, Oliver Richters. Helmholtz decomposition and potential functions for n-dimensional analytic vector fields. [10.1016/j.jmaa.2023.127138](https://doi.org/10.1016/j.jmaa.2023.127138)
