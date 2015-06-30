# hess
Code for working on conjectures related to chromatic symmetric functions and Hessenberg varieties

If you have a suitable version of GNU make installed, just run:
```
make
```

and the output should appear in the file `output.py`, a Sagemath script.
If you want things to go faster, you can do parallel execution with something like:
```
make -j maxjobs -l maxload
```

If you want to compute for a different set of Dyck path sizes, run something like:
```
make SIZES='1 2 3 4 5'
```

You can also run tests or clean the directory of all intermediate and output files:
```
make test SIZES=''
make clean SIZES=''
```
