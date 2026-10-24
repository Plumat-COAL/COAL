# V-3

Steps (Linux, either engine):

    coal big.pep big.nuc > out.fa &
    kill -9 $!
    ls $TMPDIR/.coal_*.tmp   # a leftover; SIGTERM / SIGINT / SIGHUP / a closed pipe leave none
