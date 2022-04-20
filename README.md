Small fluid flow simulation prototype in C++

* Used for experimenting fluid flow simulation in C++

This is a toy project with the building blocks for a basic Eulerian flow simulation: 

[advection] -> [add forces] -> [pressure-solver to get a divergence-free flow] -> [advect tracer particles]

I use this as a learning sandbox.

<p align="center"><img src="https://github.com/skal65535/flow_sim/blob/main/example.flow.png"></p>

The sim.cc main simulation uses SDL1 for display.
There's also a small 'gauss.cc' tool to generate the pre-iterated pressure-solver
steps code. You can change the number of steps for experimention purpose.

Code is released under MIT license.

* ShaderToy version: https://www.shadertoy.com/view/ft2czK

* Some literature:

  https://gist.github.com/vassvik/f06a453c18eae03a9ad4dc8cc011d2dc

  https://www.youtube.com/watch?v=gJMBEvYEfJI

  http://jamie-wong.com/2016/08/05/webgl-fluid-simulation/
