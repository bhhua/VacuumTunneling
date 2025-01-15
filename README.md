# VacuumTunneling

### A package for evaluating bubble profile with a renomalization factor in Mathematica

This package is designed to evaluate such an action

$$S_E=4\pi\int_{0}^{\infty} dr r^{D-1} \left[\frac{Z_{\sigma}^{-1}}{2} \left(\frac{d\sigma}{dr}\right)^2 +V_{\mathrm{eff}}(\sigma)\right],$$

whose corrsponding equation of motion reads

$$\frac{d^2 \sigma}{dr^2}+\frac{D-1}{r} \frac{d\sigma}{dr}-\frac{1}{2}\frac{\partial \mathrm{log} Z_{\sigma}}{\partial \sigma} \left(\frac{d\sigma}{dr}\right)^2 = Z_{\sigma}\frac{\partial V_{eff}}{\partial \sigma}.$$

For more details, see [蓝色的字](https://网址 "VacuumTunneling: A package to solve bonce equation with renormalization factor")

#### Installation and Running Guide

You can install this package like this

<img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/installpic1.png" width="225px"> <img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/installpic2.png" width="300px">

Load this package by

```
<<"VacuumTunneling`"
```

#### An example for single-field tunneling

Input the effective potential and renormalization factor and find the true and false vacua:

```
V[x_] := - x^2 + x^3 + 2 x^4;
Z[x_] := x^2 + 1;
tv = x /. NSolve[V'[x] == 0, x][[1]];
fv = x /. NSolve[V'[x] == 0, x][[3]];
```

Call the `Tunneling` function, input at least 5 arguments `[expression of potential, expression of renormalization, field name,
true vacuum, false vacuum]`:

```
a = Tunneling[V[x], Z[x], x, tv, fv]
```

An array of 2 elements will be outputted as the result:

<img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/basicexampleout.png" width="450px">

To extract the content, use the code

```
Plot[a[[1]][r], {r, 0, 7.57}]
Se = a[[2]]
```

<img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/basicexampleplot.png" width="300px">

#### An example for 2-field tunneling
Input the 2-field effective potential and renormalization factor and find the true and false vacua:

```
V2[h_, s_] := -0.1 h^4 -100 h^2 + 0.3 s^4 -60 s^2 +3 h^2 s^2;
Z2[h_, s_] := 0.2 / (h + 1) + 0.1 / (2 + 1) - 0.15 / ((h + 1) (s + 1));
fv2 = {x, y} /. Last[FindMinimum[V2[x, y], {{x, 0}, {y, 10}}]]
tv2 = {x, y} /. Last[FindMinimum[V2[x, y], {{x, 22.4}, {y, 0}}]]
```

Call the `Tunneling` function, note that the field names and the true/false vacua each contain two elements. You can input 1 in the place where Z is
supposed to be inputted to evaluate the tunneling without a renormalization factor.

```
b1 = Tunneling[V2[x, y], 1, {x, y}, tv2, fv2]
b2 = Tunneling[V2[x, y], Z2[x, y], {x, y}, tv2, fv2]
```

Show the bounce of each field by `Plot`:

```
Plot[{b2[[1]][x][[1]], b2[[1]][x][[2]]}, {x, 0, 2.24}]
```

<img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/2doriginplot.png" width="300px">

Show the path by `ParametricPlot`:

```
Show[
ContourPlot[v[x, y], {x, -1, 24}, {y, -1, 12}, Contours -> 50, ContourShading -> None, Epilog -> {Red, PointSize[Large], Point[{tv2, fv2}]}],
ParametricPlot[b2[[1]][x], {x, 0, 2.24}],
ParametricPlot[b1[[1]][x], {x, 0, 0.631}, PlotStyle -> {Orange, Dashed}]
]
```

<img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/2doriginpath.png" width="300px">
