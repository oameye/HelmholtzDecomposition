# HelmholtzDecomposition

Analytical Helmholtz decomposition of n-dimensional vector fields, following Richters & Glötzl, *Analytical Helmholtz Decomposition of n-Dimensional Vector Fields* ([10.1016/j.jmaa.2023.127138](https://doi.org/10.1016/j.jmaa.2023.127138)). Given a vector field `f`, the decomposition

```
f = g + r,    Curl[g] == 0,    Div[r] == 0
```

is constructed once and exposed through a queryable result object. A companion helper `HelmholtzGaugeShift` applies explicit harmonic gauge shifts to an existing decomposition.

## Install

```mathematica
PacletInstall["https://github.com/oameye/HelmholtzDecomposition/releases/download/v2.1.0/HelmholtzDecomposition-2.1.0.paclet"]
```

Or, working from a checkout:

```mathematica
PacletDirectoryLoad["/path/to/HelmholtzDecomposition"];
Needs["HelmholtzDecomposition`"];
```

## Quick start

```mathematica
Needs["HelmholtzDecomposition`"];

hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];

hd["Gradient"]                 (* {x^2 + y^2/2, x y} *)
hd["Rotational"]               (* {-(1/2)(-2 + y) y, 0} *)
hd["Potential"]                (* scalar P with Grad[P] == g *)
hd["Verify"]                   (* <|GradientCurlFree -> True, ...|> *)

List @@ hd                     (* {g, r}, handy for destructuring *)

shifted = HelmholtzGaugeShift[hd, x + y];
shifted["Gradient"]            (* g + {1, 1} *)
shifted["CanonicalAssociation"]["Gradient"]   (* original canonical g *)
```

## API — one function, queryable result

```mathematica
hd = HelmholtzDecomposition[f, xvec, opts]   (* explicit coords *)
hd = HelmholtzDecomposition[f, opts]         (* auto-detect coords from f *)
```

The result has head `HelmholtzDecomposition` and renders as a collapsible summary box in notebooks. It is queried by string property:

| Property | Meaning |
|---|---|
| `"Gradient"` | curl-free component `g` |
| `"Rotational"` | divergence-free component `r` |
| `"Potential"` | scalar `P` with `Grad[P, xvec] == g` |
| `"AntisymmetricPotential"` | the n-dimensional 2-form `Rij` whose row-divergence is `r` |
| `"VectorPotentialA"` | (3D only) vector `A` with `Curl[A] == r` |
| `"PotentialMatrix"` | the underlying potential matrix `Fij` |
| `"GaugeShift"` | harmonic vector field shifted from `r` into `g` |
| `"GaugePotential"` | scalar potential whose gradient equals `"GaugeShift"` |
| `"CanonicalAssociation"` | ungauged underlying Association |
| `"Residual"` | `f - g - r`, simplified |
| `"Verify"` | `<|"GradientCurlFree", "RotationalDivergenceFree", "ResidualZero"|>` |
| `"Field"`, `"Coordinates"` | the inputs |
| `"Properties"` | the list of valid property names |
| `"Association"` | the raw underlying Association |

`List @@ hd` returns `{g, r}` (handy for destructuring), `Normal[hd]` returns the underlying Association, `Length[hd]` returns the field dimension. Querying an unknown property returns `$Failed` with a message.

### Options

```mathematica
HelmholtzDecomposition[f, xvec,
  SimplificationFunction -> FullSimplify,   (* or Simplify, Identity, ... *)
  Assumptions :> $Assumptions,              (* forwarded to the simplifier *)
  Verbose -> False,                         (* trace the algorithm via Echo *)
  Caching -> True,                          (* memoize on (f, xvec, simp, asm, tc) *)
  TimeConstraint -> Infinity                (* circuit breaker forwarded to the simplifier *)
]
```

`Identity` is special-cased so `SimplificationFunction -> Identity` skips simplification cleanly.

## Gauge shifts

```mathematica
shifted = HelmholtzGaugeShift[hd, phi]   (* phi harmonic scalar *)
shifted = HelmholtzGaugeShift[hd, h]     (* h harmonic vector field *)
```

This returns a new `HelmholtzDecomposition` object with

```mathematica
g' = g + h
r' = r - h
P' = P + phi
```

where `h = Grad[phi, xvec]`. The helper verifies the gauge input is harmonic:

- scalar form: `Laplacian[phi, xvec] == 0`
- vector form: `Curl[h] == 0` and `Div[h] == 0`

The original canonical decomposition is preserved under
`shifted["CanonicalAssociation"]`.

Note: after a gauge shift, the matrix-level properties `"PotentialMatrix"` and `"AntisymmetricPotential"` are marked unavailable, because the package does not currently construct a closed-form coexact potential for the shifted rotational field.

### Caching

Results are memoized by default on the tuple `(f, xvec, SimplificationFunction, Assumptions)`. Repeated calls cost effectively nothing. Drop the cache with `ClearHelmholtzCache[]` (returns the number of evicted entries), or pass `Caching -> False` for a single bypass.

### Failure modes

The function returns `$Failed` (with a typed message under `HelmholtzDecomposition::*`) on:

- length mismatch between `f` and `xvec` (`::nomatch`)
- ambiguous coordinate auto-detection (`::autocoord`)
- atomic terms outside the polynomial / exponential / sine / cosine grammar (`::nodec`)
- query for an unknown property (`::nokey`)
- query for `"VectorPotentialA"` outside 3D (`::baddim`)
- non-boolean `Verbose` (`::badverbose`)
- invalid harmonic gauge input to `HelmholtzGaugeShift` (`HelmholtzGaugeShift::badharm`, `::badlen`, `::badpot`)

## Notes

**Simplifier contract.** `SimplificationFunction` must accept the call shape `simp[expr, asm, TimeConstraint -> tc]`. `FullSimplify`, `Simplify`, and the special-cased `Identity` all satisfy this. A custom simplifier `myFn[expr_] := ...` would receive arguments it doesn't expect; wrap with `Function[{e, _, _}, myFn[e]]` if needed.

**Cache key is the symbolic expression.** The cache is keyed on `(f, xvec, simp, asm, tc)` as written, not on currently-bound values. If you call `HelmholtzDecomposition[{x^2}, {x}]`, then assign `x = 5`, then call again, the cached object is returned — even though `f` would now evaluate to `{25}`. This is the right semantic (the field is a symbolic expression, parameters bind at access) but worth knowing.

**Property access re-evaluates.** `hd["Field"]` returns the stored expression and re-evaluates against current bindings. So `hd["Field"]` after `x = 5` evaluates to a numeric vector. Same caveat.

## Citation and license

Please cite the original paper and Zenodo repository when using this code.

- [Zenodo repository](https://zenodo.org/records/7680297)
- Erhard Glötzl, Oliver Richters. *Helmholtz decomposition and potential functions for n-dimensional analytic vector fields*. [10.1016/j.jmaa.2023.127138](https://doi.org/10.1016/j.jmaa.2023.127138)
