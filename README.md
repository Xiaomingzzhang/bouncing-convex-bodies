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
<video width="320" height="240" controls>
  <source src="rotate_ellispe.mp4" type="video/mp4">
</video>

