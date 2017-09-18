
python header.py
for fileprefix in $( ls *ratios*|cut -f1 -d\ |uniq); do
  cat *$fileprefix*gcl_rpe_ratios* | python process.py $fileprefix
done

python header.py
for fileprefix in $( ls *ratios*|cut -f1 -d\ |uniq); do
  cat *$fileprefix*gcl_values* | python process.py gcl_values$fileprefix
done

python header.py
for fileprefix in $( ls *ratios*|cut -f1 -d\ |uniq); do
  cat *$fileprefix*rpe_values* | python process.py rpe_values$fileprefix
done
