# Theoretical protocol

The idea of this folder is to show you how we analyse the acoustic physics of a flute. First of all, we will start by taking a few things for granted. The first one ins known as Pressure-gradient force. The second is the fact that a velocity gradient produces a compression in the fluid. Of course, this things as we shown in (1)  has to be demostrated. We will explain this with more detail in the extensive work. So as not to complicate things, we start in (1)

<img src="https://render.githubusercontent.com/render/math?math= \frac{\partial p}{\partial x} = -\rho \frac{\partial u}{\partial t} \;\;\;\;\; ; \;\;\;\;\; \kappa \frac{\partial p}{\partial t} = -\frac{\partial u}{\partial x}"> (1)

This lead to a wave equation like (2).

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/onda.png) (2)

The solution for this equation (1) it is usually as follows:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/solonda.png) (3)

Using (1) and replacing _*p*_ in (2) we get the fluid velocity:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/u.png) (4)

Multiplying (4) by the cross section _*S*_ we get the acoustic flux:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/flujo.png) (5)

However, this is for infinite dimensions. If we want a more realistic model is needed to add a reflected wave from an open or closed edge.

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/pr.png) (6)

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/ur.png) (7)

With the flux we can define the acoustic impedance in terms of _*ρ*_, _*c*_ y _*S*_ :

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/impe.png) (8)

We can figure out a lot of interesting acoustic properties of the flute. Most scientific studies base their analysis on the study of impedance. The problem with doing this kind of analisys in that the impedance measurement without a laboratory is really complex. In case of needed, we did some graphics in the folder avobe called Impedance. At the end we will tell how we did those graphics. In the meantime we will do the analisys in terms of the intensity defined by:


![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/I.png) (9)

Being _*Z'=ZS*_. Knowing that _*Z=P/U_* we have:

![](https://github.com/saguileran/Acoustics-Instruments/blob/master/Theory/Equations/Iz.png) (10)

<img src="https://render.githubusercontent.com/render/math?math= e^{i \pi} = -1">
