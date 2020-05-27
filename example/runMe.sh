#!/bin/bash

a="log.txt"
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

function compare_log_files () {

    echo "[TEST] Comparing $1 $2."

    if diff <( sort $1 ) <( sort $2 ) > /dev/null 2>&1
    then
        printf "${GREEN}[PASS]${NC} $3.\n"
    else
        printf "${RED}[FAIL]${NC} $3.\n"
    fi
    rm log.txt outlier_1.txt outlier_2.txt -f
}

rm log.txt outlier_1.txt outlier_2.txt -f

# default mode
echo "[TEST] Test 1. Default mode."
echo "[TEST] ../checker_complete.2.pl -a EOG090X002Z.aa.fas"
../checker_complete.2.pl -a EOG090X002Z.aa.fas &>/dev/null
compare_log_files $a "results/log.txt" "Test 1"


# default w/ references
echo "[TEST] Test 2. Default mode using references."
echo "[TEST] ../checker_complete.2.pl -a EOG090X002Z.aa.fas -r EOG090X002Z.aa.refs"
../checker_complete.2.pl -a EOG090X002Z.aa.fas -r EOG090X002Z.aa.refs &>/dev/null
compare_log_files $a "results/log.2.txt" "Test 2"

# default w/ references & subjects
echo "[TEST] Test 3. Default mode using references and subjects."
echo "[TEST] ../checker_complete.2.pl -a EOG090X002Z.aa.fas -r EOG090X002Z.aa.refs -s EOG090X002Z.aa.subj"
../checker_complete.2.pl -a EOG090X002Z.aa.fas -r EOG090X002Z.aa.refs -s EOG090X002Z.aa.subj &>/dev/null
compare_log_files $a "results/log.3.txt" "Test 3"

