
python header.py
for fileprefix in $( ls *ratios*|cut -f1 -d\ |uniq); do
  cat *$fileprefix*gcl_rpe_ratios* | python process.py $fileprefix
  done
