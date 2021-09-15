#!/bin/bash

sbatch -p general -N 1 -n 16 -t 60:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r run_analysis\(\) -logfile r01.log"

