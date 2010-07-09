#!/bin/bash

pep=ext_and_int_altspliced.sorted.nopseudo.pep
# LS_SS_PCL.diff_internal.ids.nopseduo.pep

pepshortid=${pep%%pep}shortid.pep
perl -ne 'chomp; if(m/^>(\S+)/) { print ">$1\n"; } else { print $_,"\n"; }' < $pep > $pepshortid
../src/phobius/phobius.pl -long ${pep%%pep}shortid.pep > ${pepshortid}.phobius_long
../nematodes/bin/mask_phobius_signal.pl ${pepshortid} ${pepshortid}.phobius_long> ${pep%%pep}phobius_signal_masked.pep

