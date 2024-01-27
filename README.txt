-----------------------------------------------------------
 How to use the derivative-free optimizer SDMINMAX
-----------------------------------------------------------

1- Gunzip and untar the archive in a folder on your computer.

2- Edit file problem.f to define your own problem.
   In particular, you need to define the routines
   set_dimension 	which sets the dimensions of the problem
   startp	 	which sets the starting point
   functs 		which returns the component function values

3- At command prompt execute 

     $> make
 
   which will create the executable 'minmax'

4- execute

     $> ./minmax
