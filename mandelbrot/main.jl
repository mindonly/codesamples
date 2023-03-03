#!/usr/bin/env julia

include("mandelbrot.jl")

# left side real line intersection
is_mandelbrot(-1.99)
is_mandelbrot(-2)
is_mandelbrot(-2.01)
is_mandelbrot(-2.001)
# right side real line intersection
is_mandelbrot(0.249)
is_mandelbrot(0.25)
is_mandelbrot(0.251)
# show sequence
is_mandelbrot(0.251, sq=true)
is_mandelbrot(0.2501)
is_mandelbrot(0.25001)
is_mandelbrot(0.250001)
# false positive
is_mandelbrot(0.2500001, it=4000)
is_mandelbrot(0.2500001, it=10000)
# geogebra webapp
is_mandelbrot(-1 + 0.2864im)
is_mandelbrot(-1 + 0.2865im)
is_mandelbrot(-1 + 0.2865im, sq=true)
# another example
is_mandelbrot(-1 + 0.28646im)
is_mandelbrot(-1 + 0.28647im)

# construct Mandelbrot set
setup_mandelbrot()
setup_mandelbrot(r=2.5, px=500, cs=1, ra=π/6)
setup_mandelbrot(r=2.25, px=800, cs=1, yr=true, dr='v', ra=-π/4)

# alias plots
 left()
   up(r=3)
right()
 down(r=1)
# gradients
plot_gradients(iter=16)
plot_gradients(cs=1, iter=32)

# rotations
right()
right(ra= π/6)
right(ra= π/4)
right(ra= π/2)
right(ra=5π/6 + π, cs=0)
right(ra= π/4 - 4π/3, cs=0)
# unit circle
rotate_unit_circle(r=2, cs=1)
rotate_unit_circle(r=1)
