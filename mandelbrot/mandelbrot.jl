#!/usr/bin/env julia

using Images

# quadratic map: z, c complex
f(z, c) = z^2 + c

#  determine Mandelbrot set membership for a complex number c
#     returns boolean T/F
function is_mandelbrot(c; sq=false, it=4000)
   MAX_ITER = it 
   z = 0
   println("\nc = $c is an element of the Mandelbrot set.")
   
   for i in range(1, MAX_ITER)
      z = f(z, c)
      mag = abs(z)      # modulus, magnitude
      
      if sq == true
         @show i, z, mag
      end
      if mag > 2        # who proved this condition?
         if sq == false
            @show i, z, mag
         end
         return false
      end
   end
   
   return true
end

#  complex parameter c, max.iterations iter
#     returns ratio iter.count / max.iterations
function mandelbrot(c, iter)
   MAX_ITER = iter 
   n = 0    # iter count
  
   # iterated quadratic map
   z = 0
   while abs(z) <= 2 && n < MAX_ITER
      z = z^2 + c
      n += 1
   end
   
   return n / MAX_ITER
end

#  plot square radius r, canvas pixels px, y-reflected yr (t/f), rot. angle ra (radians)
#  color scheme cs (0 or 1), iterations it, direction dr (horiz. 'h' or vert. 'v')
#     returns Mandelbrot set image data pixel matrix
function setup_mandelbrot(; r=2, px=1000, yr=false, ra=0, cs=0, it=32, dr='h')
   # plot window, real & imaginary axes
   RE_MIN, RE_MAX = [-r, r]
   IM_MIN, IM_MAX = [-r, r]
   
   # canvas size (pixels) max.1000
   if px > 1000; px = 1000; end
   HEIGHT, WIDTH = [px, px]
   
   # initialize w/ random image data (grayscale objects)
   img_data = rand(Gray, HEIGHT, WIDTH)

   # set up pixel matrix
   for x in range(1, WIDTH)
      for y in range(1, HEIGHT)
         # reflect over y-axis
         if yr == true
            sgn = -1
         else sgn = 1
         end

         # convert pixel coords to complex
         c = complex(sgn*RE_MIN + (x / WIDTH)  * sgn*(RE_MAX - RE_MIN),
                     sgn*IM_MIN + (y / HEIGHT) * sgn*(IM_MAX - IM_MIN))
         # rotation angle
         θ = ra
         c = c * ℯ^(im*θ)

         # color scheme
         if cs == 1
            pixel = mandelbrot(c, it)       # white/black
         else
            pixel = 1 - mandelbrot(c, it)   # black/white
         end
         # image data update direction (horiz or vert)
         if dr == 'v'
            img_data[x, y] = Gray(pixel)
         else
            img_data[y, x] = Gray(pixel)
         end
      end
   end
   
   return img_data
end

# setup_mandelbrot() aliases
 left(;r=2,
      px=1000,
      cs=1,
      it=32,
      ra=0) = setup_mandelbrot(r=r, px=px, cs=cs, it=it, ra=ra)
   up(;r=2,
      px=1000,
      cs=0,
      it=32,
      ra=0) = setup_mandelbrot(r=r, px=px, cs=cs, it=it, ra=-ra, dr='v')
right(;r=2,
      px=1000,
      cs=1,
      it=32,
      ra=0) = setup_mandelbrot(r=r, px=px, cs=cs, it=it, ra=ra, yr=true)
 down(;r=2,
      px=1000,
      cs=0,
      it=32,
      ra=0) = setup_mandelbrot(r=r, px=px, cs=cs, it=it, ra=-ra, yr=true, dr='v')

#  gradients
function plot_gradients(; cs=0, iter=32)
   for j in range(1, iter)
      display(left(cs=cs, it=j))
      sleep(0.5)
   end
end

#  high school unit circle
function rotate_unit_circle(; r=2, cs=0)
   pi6 = [ m*π/6 for m in 1:12 ]
   pi4 = [ n*π/4 for n in 1:8 ]
   angles = vcat(pi6, pi4)
   push!(angles, 0)
   angles = unique(sort!(angles))
   
   for t in angles
      display(right(r=r, cs=cs, ra=t))
      sleep(0.35)
   end
end

# https://www.codingame.com/playgrounds/2358/how-to-plot-the-mandelbrot-set/mandelbrot-set
# https://realpython.com/mandelbrot-set-python/
