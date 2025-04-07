# Numerical simulations of bouncing convex bodies
Numerical simulations for the paper: Integrability and chaos in the bouncing convex body model.

## Verification of the formulas of the lowest point when rotating the ellipse and the eccentricity disk:
Mathematica code for ellipse:
```mathematica
With[{a = 2, b = 1},
 Animate[
  Graphics[{Rotate[Disk[{0, 0}, {a, b}], phi], Red, PointSize[Large], 
    Point[{((-a^2 + b^2)  Cos[phi]  Sin[phi])/Sqrt[
      b^2  Cos[phi]^2 + a^2  Sin[phi]^2], -Sqrt[
       b^2  Cos[phi]^2 + a^2  Sin[phi]^2]}]}, 
   PlotRange -> {{-3, 3}, {-3, 3}}], {phi, 0, 2 Pi}]]
```

<img src="rotate_ellipse.gif" width="500"/>


Mathematica code for eccentricity disk:

```mathematica
With[{r = 2, c = 1},
 Animate[
  Graphics[{Rotate[Circle[{c, 0}, r], phi, {0, 0}], Black, 
    PointSize[Large], Point[{0, 0}], 
    Rotate[Point[{c, 0}], phi, {0, 0}], Red, PointSize[Large], 
    Point[{c  Cos[phi], -r + c  Sin[phi]}]}, 
   PlotRange -> {{-4, 4}, {-4, 4}}], {phi, 0, 2 Pi}]]
```

<img src="rotate_circle.gif" width="500"/>

## Symbolic computation for $h_{22}(\phi_{i-1},\phi_{i})+h_{11}(\phi_{i},\phi_{i+1})$

Mathematica code:
```mathematica

```
