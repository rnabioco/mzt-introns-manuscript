#! /usr/bin/env bash

rm -f trackDb.txt 

bash dms/makeDmsDb.sh > dms/trackDb.txt
bash enzymatic_probing/makeEpDb.sh > enzymatic_probing/trackDb.txt

dbfiles=$(find . -maxdepth 2 -name 'trackDb.txt')
cat tDb.tmp $dbfiles >> trackDb.txt

