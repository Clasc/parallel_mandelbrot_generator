# OpenMP assignment

to build the project run:
```./build.sh ```

to run the project, you can pass in 2 variables. The first one is the number of threads and the second is the number of tasks per thread. If none are given the default is 16 threads and 4 tasks per thread (this has in my experience the best performance):
```./out/openmp <threads> <tasks_per_thread> ```
example:
```./out/openmp 16 4 ```.

If you are on vscode, you can launch the project by clicking F5 (debug), or ctrl F5 for no debugging.
The script ./build-alma.sh is used to build the project on alma with the corresponding compiler version.

After everything is run, you can run: 
```./run_n_times.sh <iterations> ```
to run the program n times.

This is useful to collect some data. It creates a file out.txt and writes the log to it.
If the out.txt file already exists, it will be overwritten.

When everything is run, you can run
```./assert.sh ```
to test the output.
This script checks if the diff of the created mandelbrot ppm file and the original one is empty.
OK if it is empty
NOK if it is not.

*Warning*
The test will be also "OK" if the ppm file does not exist. So please make sure that you run them at lest once!

