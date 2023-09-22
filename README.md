<h1>Introduction</h1>

<pre>
c++ -O2 main.cpp -lglut -lGL
</pre>

or
<pre>
c++ -O2 main.cpp -lglut `pkg-config --libs gl`
</pre>

Run
<pre>
./a.out 1 8 3
./a.out 1 8 3 saveimages
</pre>

<p align="center"><img src="img/cover.png"/></p>

<h1>References</h1>

<a href="https://en.wikipedia.org/wiki/Granular_convection">wikipedia:Granular_convection</a>

<h1>Method</h1>

a collision step between two particles</br>
</br>
INPUT: dt, radius, relative position r, v1, v2, omega1, omega2 are the
velocities and angular velocities of the two colliding particles</br>
</br>
OUTPUT: force1, force2, torque1, torque2

<h1>Results</h1>

<p align="center"><img src="img/box.gif"/></p>
