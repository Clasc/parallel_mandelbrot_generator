# OpenMP 

to build the project run:
```./build.sh ```

to run the project, you can pass in 2 variables. The first one is the number of threads and the second is the number of tasks per thread. If none are given the default is 16 threads and 4 tasks per thread (this has in my experience the best performance):
```./out/openmp [] ```