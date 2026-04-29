# HelmholtzDecomposition

A Mathematica package for the analytical Helmholtz decomposition of vector fields in n dimensions, following Richters & Glötzl, *Analytical Helmholtz Decomposition of n-Dimensional Vector Fields* ([10.1016/j.jmaa.2023.127138](https://doi.org/10.1016/j.jmaa.2023.127138)). Given `f`, the decomposition

```
f = g + r,    Curl[g] = 0,    Div[r] = 0
```

is constructed from a potential matrix `Fij`, with the scalar potential `P` (such that `g = ∇P`) and the antisymmetric matrix `Rij` (whose row-divergence is `r`) made directly available.

## Install

```mathematica
PacletInstall["https://github.com/oameye/HelmholtzDecomposition/releases/download/v1.1.0/HelmholtzDecomposition-1.1.0.paclet"]
```

Or, working from a checkout:

```mathematica
PacletDirectoryLoad["/path/to/HelmholtzDecomposition"];
Needs["HelmholtzDecomposition`"];
```

## Quick start

```mathematica
Needs["HelmholtzDecomposition`"];

f = {x^2 + y, x y};

(* Coordinates auto-detected from f *)
{g, r} = HelmholtzDecomposition[f];

(* Scalar potential P with Grad[P] == g *)
P = HelmholtzPotential[f];

(* Everything in one shot, computed once and cached *)
data = HelmholtzAssociation[f];
data["Potential"]              (* P *)
data["VectorPotentialMatrix"]  (* Rij *)
data["Residual"]               (* should be {0, 0} *)
```

## API

| Function | Returns |
|---|---|
| `HelmholtzDecomposition[f, xvec]` | `{g, r}` (gradient and rotational components) |
| `HelmholtzGradient[f, xvec]` | `g` (curl-free component) |
| `HelmholtzRotational[f, xvec]` | `r` (divergence-free component) |
| `HelmholtzPotential[f, xvec]` | scalar `P` with `Grad[P, xvec] == g` |
| `HelmholtzVectorPotential[f, xvec]` | antisymmetric matrix `Rij` whose row-divergence is `r` |
| `HelmholtzAssociation[f, xvec]` | `Association` with all of the above plus `"PotentialMatrix"`, `"Residual"`, `"Coordinates"`, `"Field"` |
| `PrintHelmholtz[f, xvec]` | prints a diagnostic block **and returns the same Association** |
| `PrintHelmholtz[assoc]` | prints diagnostic for an already-computed Association (no recomputation) |

The `xvec` argument is optional in every signature: omit it and the package extracts the coordinate symbols from `f`. Pass it explicitly when the field contains parameter symbols you don't want treated as coordinates.

### Options

All surface functions accept the following options:

- `Assumptions -> ...` — forwarded to the simplifier (`$Assumptions` by default). Important for parameterised fields where signs/positivity matter.
- `Simplify -> FullSimplify | Simplify | Identity` — controls the simplifier applied to `g`, `r`, `P`, `Rij`. `Identity` skips simplification (fast, but verbose output).

```mathematica
HelmholtzPotential[{a x^2, b y}, {x, y}, Assumptions -> a > 0 && b > 0]
HelmholtzAssociation[bigField, coords, Simplify -> Identity]
```

### Caching

`HelmholtzAssociation[f, xvec, opts]` is memoized on the `(f, xvec, simplifier, assumptions)` tuple, so calling `HelmholtzPotential` and `HelmholtzGradient` back-to-back on the same field does the heavy `Fij` computation only once. All surface functions go through `HelmholtzAssociation` internally.

### Verbose mode

```mathematica
SetVerbose[True];   (* trace which atoms hit which branch of AutocalculateFij *)
GetVerbose[]
```

## Citation and license

Please cite the original paper and Zenodo repository when using this code.

- [Zenodo repository](https://zenodo.org/records/7680297)
- Erhard Glötzl, Oliver Richters. *Helmholtz decomposition and potential functions for n-dimensional analytic vector fields*. [10.1016/j.jmaa.2023.127138](https://doi.org/10.1016/j.jmaa.2023.127138)
