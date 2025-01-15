# VacuumTunneling

### A package for evaluating bubble profile with a renomalization factor in Mathematica

This package is designed to evaluate such an action

$$S_E=4\pi\int_{0}^{\infty} dr r^{D-1} \left[\frac{Z_{\sigma}^{-1}}{2} \left(\frac{d\sigma}{dr}\right)^2 +V_{\mathrm{eff}}(\sigma)\right],$$

whose corrsponding equation of motion reads

$$\frac{d^2 \sigma}{dr^2}+\frac{D-1}{r} \frac{d\sigma}{dr}-\frac{1}{2}\frac{\partial \mathrm{log} Z_{\sigma}}{\partial \sigma} \left(\frac{d\sigma}{dr}\right)^2 = Z_{\sigma}\frac{\partial V_{eff}}{\partial \sigma}.$$

For more details, see [蓝色的字](https:// "VacuumTunneling: A package to solve bonce equation with renormalization factor")

#### Installation and Running Guide

You can install this package like this

<img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/installpic1.png" width="225px"> <img src="https://github.com/bhhua/VacuumTunneling/blob/main/images/installpic2.png" width="300px">

Load this package by

```
<<"VacuumTunneling`"
```

#### An example for single-field tunneling

Input the effective potential and renormalization factor:

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

```
Plot[a[[1]][r], {r, 0, 8.01}]
Se = a[[2]]
```
