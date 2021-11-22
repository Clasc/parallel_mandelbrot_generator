#usr/bin/bash
difference="$(diff mandelbrot.ppm Mandelbrot_parallel.ppm)"
[[ ! -z "$difference" ]] && echo "NOK" || echo "OK"
