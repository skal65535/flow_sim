Small C++ prototype of fluid flow simulation.

<p align="center"><img src="https://github.com/skal65535/flow_sim/blob/main/example.flow.png"></p>

This is a toy project for playing with the building blocks of a basic Eulerian flow simulation: 

[advection] -> [add forces] -> [pressure-solver to get a divergence-free flow] -> [advect tracer particles]

Code is released under MIT license.

# ShaderToy version:

   https://www.shadertoy.com/view/ft2czK

# Building

src/sim.cc is the simulation code. It can display the flow using SDL1.
There's also a small 'gauss.cc' tool to generate the pre-iterated pressure-solver
steps code, in case you want to change the number of Jacobi steps in the
code.

You can build with 'cmake':

   `mkdir build && cd build && cmake ../ && make -j`

# Some literature:

  https://gist.github.com/vassvik/f06a453c18eae03a9ad4dc8cc011d2dc

  https://www.youtube.com/watch?v=gJMBEvYEfJI

  http://jamie-wong.com/2016/08/05/webgl-fluid-simulation/
