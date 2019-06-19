#!/bin/bash

#inits
export ASPERA_SCP_PASS="z1BK7fEqO5OL";
username="dbox2054"
destination="`dirname $0`"
parallel_downloads=8

#overwrite default parameters?
if [ ! -z $1 ]; then parallel_downloads="$1"; fi
if [ ! -z $2 ]; then destination="$2"; fi


#get small files in its most convenient way
/home/Shared_taupo/steph/src/.aspera/cli/bin/ascp --ignore-host-key -E "*.aes" -E "*.cip" -E "*.crypt" -E "download.sh" -d -QTl 100m ${username}@xfer.crg.eu: ${destination}/

#get not small files in its most convenient way
cat ${destination}/dbox_content_file.txt | xargs -i --max-procs=$parallel_downloads bash -c "mkdir -p $destination/\`dirname {}\`; echo \"Downloading {}, please wait ...\"; /home/Shared_taupo/steph/src/.aspera/cli/bin/ascp --ignore-host-key -k 1 --partial-file-suffix=PART -QTl 100m ${username}@xfer.crg.eu:{} ${destination}/{} >/dev/null 2>&1"

