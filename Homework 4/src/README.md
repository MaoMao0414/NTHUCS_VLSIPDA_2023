# Homework 4 Global Placement
Follow the tips in `GlobalPlacer.cpp` to get a placement with minimum HPWL.

## How to Compile
In this directory, enter the following command:
```sh
$ make
```
It will generate the executable file `hw4` in the directory `HW4/bin`.
And it will also generate a executable file of parellal version named as `hw4_parallel` in the directory `HW4/bin`.

If you want to remove it, please enter the following command:
```sh
$ make clean
```

## How to Run
In bin directory
Usage:
```
$ ./hw4 <.aux file> <.gp.pl file>
$ ./hw4_parallel <.aux file> <.gp.pl file>
```
For example,
```sh
$ ./hw4 ../testcase/public1/public1.aux ../output/public1.gp.pl
$ ./hw4_parallel ../testcase/public1/public1.aux ../output/public1.gp.pl
```

## How to Test
## Note: I don't implement this function for parellal version
In this directory, enter the following command:
```sh
$ make test ${name}
```
It will test on ${name} and verify the result.

For example, test on public1.
```sh
$ make test public1
```
