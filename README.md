# Numerical simulations of bouncing convex bodies

## Verification of the formulas of the lowest point when rotating the ellipse and the eccentricity circle:
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

<video src="rotate_ellipse.mp4" width="450"  controls></video>

Mathematica code for eccentricity circle:

```mathematica
With[{r = 2, c = 1},
 Animate[
  Graphics[{Rotate[Circle[{c, 0}, r], phi, {0, 0}], Black, 
    PointSize[Large], Point[{0, 0}], 
    Rotate[Point[{c, 0}], phi, {0, 0}], Red, PointSize[Large], 
    Point[{c  Cos[phi], -r + c  Sin[phi]}]}, 
   PlotRange -> {{-4, 4}, {-4, 4}}], {phi, 0, 2 Pi}]]
```

<video src="rotate_circle.mp4" width="450"  controls></video>