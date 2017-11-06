# SCAHA: Structure-based Chemical Shift Assignments via the Hungarian Algorithm
  
- SCAHA is implemented as an R script

## Install Dependencies (in an R shell)
```shell
source("R/library.R")
install_scaha_dependencies()
```

## Usage Info
```shell
$ bin/scaha -h

Usage: ./bin/scaha [options] predicted_chemical_shifts chemical_shift_peaks

Options:
	-p, --parallel
		run assignment in parallel [default FALSE]

	-n NPROCESSORS, --nprocessors=NPROCESSORS
		number of processor to use for parallel processing [default 18]

	-o OUTPUT, --output=OUTPUT
		output prefix name [default tests/assigned_shifts]

	-t, --testing
		run in testing mode; correct assignment are known and provided [default FALSE]

	-v, --verbose
		print header and progress information [default FALSE]

	-h, --help
		Show this help message and exit
```

## Input (computed chemical shifts file)
### format
_conformation, residue-number, residue-name, nucleus-type, cs-value, id-tag_
### example
```shell
1 1 GUA C5' 65.9103 2KOC
1 1 GUA C4' 83.1202 2KOC
1 1 GUA C3' 73.7023 2KOC
1 1 GUA C2' 75.5832 2KOC
1 1 GUA C1' 92.2671 2KOC
1 1 GUA C8 138.81 2KOC
1 1 GUA N1 146.484 2KOC
1 1 GUA H5' 4.42876 2KOC
1 1 GUA H5'' 4.26016 2KOC
1 1 GUA H4' 4.51489 2KOC
1 1 GUA H3' 4.62198 2KOC
1 1 GUA H2' 4.70497 2KOC
1 1 GUA H1' 5.7441 2KOC
1 1 GUA H8 8.01219 2KOC
1 1 GUA H1 12.5966 2KOC
...
...
```

## Input (assigned observed chemical shifts file)
### format
_residue-name, residue-number, nucleus-type, cs-value, dummy-value_
### example
```shell
GUA 1 H1 13.01 .
GUA 1 H1' 5.86 .
GUA 1 H2' 4.98 .
GUA 1 H3' 4.46 .
GUA 1 H4' 4.6 .
GUA 1 H5' 4.43 .
GUA 1 H5'' 4.33 .
GUA 1 H8 8.17 .
GUA 1 C1' 90.9 .
GUA 1 C2 156.3 .
GUA 1 C2' 75 .
GUA 1 C3' 74.2 .
GUA 1 C4 153.1 .
GUA 1 C4' 83.7 .
GUA 1 C5 118.8 .
GUA 1 C5' 67.7 .
GUA 1 C6 157.1 .
GUA 1 C8 139 .
GUA 1 N1 148.6 .
GUA 1 N3 161.7 .
GUA 1 N7 233.1 .
GUA 1 N9 169.5 .
GUA 1 P -1.09 .
...
...
```

## Input (truly unassigned observed chemical shifts file)
### format
_cs-value_
### example
```shell
13.01
5.86
4.98
4.46
4.6
4.43
4.33
8.17
90.9
156.3
75
74.2
153.1
83.7
118.8
67.7
157.1
139
148.6
161.7
233.1
169.5
-1.09
...
...
```

## Examples
```shell
$ # assign chemical shift contained in an assigned chemical shift data file based on LARMORD computed chemical shifts 
$ # in testing mode the actual assigned chemical shift is included in the output
./bin/scaha tests/larmord_2KOC_single.txt tests/observed_shifts_assigned_2KOC.txt --verbose --testing --output=tests/assigned_shifts
$
$ # assign chemical shift contained in a simple list of unassigned chemical shift peaks
$ # in this mode, tests/observed_shifts_assigned_2KOC.txt should contain a single column of observed peak values
./bin/scaha tests/larmord_2KOC_single.txt tests/observed_shifts_unassigned_peaks_2KOC.txt --verbose --output=tests/assigned_shifts --output=tests/assigned_shifts
$
$ # assign chemical shift contained in an assigned chemical shift data file based on LARMORD chemical shifts computed from a set of conformations (here 30 conformations)
./bin/scaha tests/larmord_2KOC_pool.txt tests/observed_shifts_assigned_2KOC.txt --verbose --output=tests/assigned_shifts
$
$ # same as above, except SCAHA is being executed in parallel mode, 
$ # which is efficient if being executed in a mult-thread computing environment 
$ # and useful when assigning chemical shifts to more than 1 conformation 
./bin/scaha tests/larmord_2KOC_pool.txt tests/observed_shifts_assigned_2KOC.txt --verbose --output=tests/assigned_shifts --parallel --nprocessors=4
```

