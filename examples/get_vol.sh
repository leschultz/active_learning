#!/bin/bash
cat $1 | tail -n 1 | awk '{print $2}' > new_vol.txt
echo "$(echo -n 'variable newVol equal '; cat new_vol.txt)" > new_vol.txt
