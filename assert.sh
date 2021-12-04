#usr/bin/bash
echo "Asserting Tasks solution..."
difference="$(diff src/lib/a2/mandelbrot-original.ppm Mandelbrot_parallel.ppm)"
[[ ! -z "$difference" ]] && echo "NOK" || echo "OK"

echo "Asserting Taskloops solution..."
difference="$(diff src/lib/a2/mandelbrot-original.ppm Mandelbrot-taskloop.ppm)"
[[ ! -z "$difference" ]] && echo "NOK" || echo "OK"
