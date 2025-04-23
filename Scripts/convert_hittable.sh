#!/bin/bash

cut -f2 westnile_hit.txt | sort | uniq > westnile_accessions.txt
grep -v '^#' westnile_accessions.txt | grep -v '^$' > westnile_accessions_cleaned.txt

