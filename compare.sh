#usr/bin/bash
difference="$(diff src/lib/a2/mandelbrot-original.ppm Mandelbrot_parallel.ppm)"
[[ ! -z "$difference" ]] && echo "NOK" || echo "OK"
