#!/bin/bash
while true; do
    date >> ./dates.log
    find /scratch/twlab/yliang/ -exec touch {} +
    sleep 2678400
done

## find will execute grep and will substitute {} with the filename(s) found. 
# The difference between ; and + is that with ; a single grep command for each file is executed 
# whereas with + as many files as possible are given as parameters to grep at once.