# Theoretical protocol

The idea of this folder is to show you how we analyse the acoustic physics of a flute. First of all, we will start by taking a few things for granted. The first one is known as Pressure-gradient force. The second is the fact that a velocity gradient produces a compression in the fluid. Of course, this things as we show in (1)  have to be demostrated. We will explain this with more detail in the extensive work. So as not to complicate things, we start in (1)

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/gradiente.png) (1)

This lead to a wave equation like (2).

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/onda.png) (2)

The solution for this equation (1) is usually as follows:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/solonda.png) (3)

Using (1) and replacing _*p*_ in (2) we get the fluid velocity:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/u.png) (4)

Multiplying (4) by the cross section _*S*_ we get the acoustic flux:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/flujo.png) (5)

However, this is for infinite dimensions. If we want a more realistic model, it is necessary to add a reflected wave from an open or closed edge.

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/pr.png) (6)

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/ur.png) (7)

With the flux we can define the acoustic impedance in terms of _*ρ*_, _*c*_ y _*S*_ :

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/impe.png) (8)

We can figure out a lot of interesting acoustic properties of the flute. Most scientific studies base their analysis on the study of impedance. The problem with doing this kind of analisys is that the impedance measurement without a laboratory and specialized equipment is really complex. In needed, we uploaded some plots in the folder above called Impedance. At the end, we will explain how those plots were made. In the meantime, we will do the analisys in terms of the intensity defined by:


![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/I.png) (9)

Being _*Z'=ZS*_. Knowing that _*Z=P/U_* we have:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/ImpAbs.PNG) (10)

### Calculating U²

The most cumbersome part of the theoretical study has been the calculation of U². We start by squaring equation (7), considering a complete reflection so that:


![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/RazonAB.PNG) (11)

Meaning that


![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/u2.png) (12)

Doing some calculations and labeling the complex conjugate by (*)


![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/cc.png) (13)

Multiplying:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/ll.png) (14)

And using Euler's formula for cosine:


![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/FlujoCuadrado.PNG) (15)

Because we are considering a complete reflection, other authors have found that it is possible to do the math at x=0 and it should be the same as doing it at any point between 0<xL [1, 2]:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/f.png) (16)

Here we use the fact that k=2πf/c

### Sound intensity

Substituing (16) in (10) gives us:


![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/ImpZ.PNG) (17)

The calculations leading to the following impedance equation are quite extensive and will not be developed here. Instead, we redirect you to [1]

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/zt.png) (18)

Where *α* and *v* are adjustable parameters that depend on temperature and the flute's geometry. Given that this is a complex number we need to do:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/lzl.png) (19)

 And finally we have:
 
 ![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/ifinal.png) (20)

But most equipments measure in SIL, so we need to convert (20) from W/m² to dB units:

 ![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/db.png) (21)
 
 Here Io is a reference sound intensity wich value is 1W/m² and I is the whole equation (20).

## Plots

First, we plot the impedance from equation  (19) in a frequency range between 1 and 5000 with the code below. The values of *α* and *v* were also taken from [1]:

 ![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/av.png) (22)

```
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstdio>

double zn(double b, double d);
double zc(double g, double h);

int main()
{
const float fmin = 0.0; //Frecuencia Minima
const float fmax = 5000.0; //Frecuencia Maxima
float f = fmin; //Frecuencia
double al; // Parametro alpha 1/m
const float L = 0.2; //Longitud Flauta m
double v; //Parametro v m/s
const float c = 343.0;//Velocidad del sonido m/s
const float a = 0.02; //Radio de la flauta m
double ar1; //Argumento 1 de la funcion
double ar2; //Argumento 2 de la funcion
const float df = 10.0;
const float N = (fmax-fmin)/df;

for(int n = 1; n <= N; n++){ //Pasos de avance de frecuencia
    f = fmin + df*n;
    al = (0.00003*(std::sqrt(f)))/a;
    v = c*(1.0 - ((0.00165)/(a*(std::sqrt(f)))));
    ar1 = al*L;
    ar2 = ((f*L*2*(M_PI))/(v));
    std::printf("%f \t %25.15e\n", f, (zn(ar1, ar2))/(zc(ar1, ar2)));     
}
return 0;
}

double zn(double b, double d){
 return std::sqrt(((std::tanh(b))*(std::tanh(b)))+((std::tan(d))*(std::tan(d))));
}

double zc(double g, double h){
 return std::sqrt(1.0+((std::tanh(g))*(std::tanh(g)))*((std::tan(h))*(std::tan(h))));
}
````
This returns some data that we can plot:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Impedance/ImpedanciaParamV1.PNG) (1)

Using the impedance obtained from the code above we made another script that computes the sound intensity in a frequency range between 1 and 5000:

````
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstdio>

double zn(double b, double d);
double zc(double g, double h);
double U(double j); //Norma de U cuadrado
int main()
{
const float fmin = 0.0; //Frecuencia Minima
const float fmax = 5000.0; //Frecuencia Maxima
float f; //Frecuencia
double al; // Parametro alpha 1/m
const float L = 0.2; //Longitud Flauta m
double v; //Parametro v m/s
const float c = 343.0;//Velocidad del sonido m/s
const float a = 0.02; //Radio de la flauta m
double ar1; //Argumento 1 de la funcion
double ar2; //Argumento 2 de la funcion
double ar3; //Argumento de la funcion U
double Z; //zn/zc
const float df = 10.0; //Paso entre cada f.
const float N = (fmax-fmin)/df; //Numero de cuentas

for(int n = 1; n <= N; n++){ //Pasos de avance de frecuencia
    f = fmin + df*n;
    al = (0.00003*(std::sqrt(f)))/a;
    v = c*(1.0 - ((0.00165)/(a*(std::sqrt(f)))));
    ar1 = al*L;
    ar2 = ((f*L*2.0*(M_PI))/(v));
    ar3 = (4.0*(M_PI)*f*L)/c;
    Z = (zn(ar1, ar2))/(zc(ar1, ar2));
    std::printf("%f \t %25.15e\n", f, U(ar3)*Z);     
}
return 0;
}

double zn(double b, double d){
 return std::sqrt(((std::tanh(b))*(std::tanh(b)))+((std::tan(d))*(std::tan(d))));
}

double zc(double g, double h){
 return std::sqrt(1.0+((std::tanh(g))*(std::tanh(g)))*((std::tan(h))*(std::tan(h))));
}

double U(double j){
 return std::sqrt(6.0 + 8.0*(std::cos(j)) + 2.0*(std::cos(2.0*j)));
}
````

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Intensity/IntensidadC.png) (2)

Wich gives figure (2). The value of A was obtained by adjusting our data to the ones found by the experiments. Because we want to predict more than one note, we would have to model the finger holes. Due to boundary conditions that required very complex calculations, only the effective length of the flute was changed, making this a more effective aproach rather than one based on first principles.

Another thing that was taken into account were hidden harmonics. Because the sound produced by a flute with a certain fingerhole configuration is not a completely pure, we add the harmonics of D to C (With L=15 and 16 cm, respectively)

````
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstdio>

double zn(double b, double d);
double zc(double g, double h);
double U(double j); //Norma de U cuadrado
int main()
{
const float fmin = 0.0; //Frecuencia Minima
const float fmax = 5000.0; //Frecuencia Maxima
float f; //Frecuencia
double al; // Parametro alpha 1/m
float L1 = 0.16; //Longitud Flauta m
float L2 = 0.15; //Longitud Flauta m
double v; //Parametro v m/s
const float c = 343.0;//Velocidad del sonido m/s
const float a = 0.02; //Radio de la flauta m
double ar1L1; //Argumento 1 de la funcion L=16
double ar2L1; //Argumento 2 de la funcion L=16
double ar3L1; //Argumento de la funcion U L=16
double ar1L2; //Argumento 1 de la funcion L=15
double ar2L2; //Argumento 2 de la funcion L=15
double ar3L2; //Argumento de la funcion U L=15
double ZL1; //zn/zc
double ZL2; //zn/zc
const float df = 10.0; //Paso entre cada f.
const float N = (fmax-fmin)/df; //Numero de cuentas

for(int n = 1; n <= N; n++){ //Pasos de avance de frecuencia
    f = fmin + df*n;
    al = (0.00003*(std::sqrt(f)))/a;
    v = c*(1.0 - ((0.00165)/(a*(std::sqrt(f)))));
    ar1L1 = al*L1;
    ar2L1 = ((f*L1*2.0*(M_PI))/(v));
    ar3L1 = (4.0*(M_PI)*f*L1)/c;
    ar1L2 = al*L2;
    ar2L2 = ((f*L2*2.0*(M_PI))/(v));
    ar3L2 = (4.0*(M_PI)*f*L2)/c;
    ZL1 = (zn(ar1L1, ar2L1))/(zc(ar1L1, ar2L1));
    ZL2 = (zn(ar1L2, ar2L2))/(zc(ar1L2, ar2L2));
    std::printf("%f \t %25.15e\n", f, (U(ar3L1)*ZL1)+(U(ar3L2)*ZL2));     
}
return 0;
}

double zn(double b, double d){
 return std::sqrt(((std::tanh(b))*(std::tanh(b)))+((std::tan(d))*(std::tan(d))));
}

double zc(double g, double h){
 return std::sqrt(1.0+((std::tanh(g))*(std::tanh(g)))*((std::tan(h))*(std::tan(h))));
}

double U(double j){
 return std::sqrt(6.0 + 8.0*(std::cos(j)) + 2.0*(std::cos(2.0*j)));
}

````
This returns (3)

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Intensity/IntensidadL15yL16.png) (3)

[1]N. H. Fletcher,The Physics of Musical Instruments, 2.aed. Springer, 2008.
[2]L. E. Kinsler,Fundamentals Of Acoustics, 1.aed. Wiley ed, 1950.




