rm -rfv output.namelist*
rm -rfv run.stat*
rm -rfv ocean.output
rm -rfv NABE_*_restart*.nc
toto=`awk '/output.path = /{print $NF}' output.conf`
rm -rfv $toto
rm -rfv output_pisces/*
rm -rfv NABE_*grid_*.nc
rm -rfv NABE_*ptrc_T*.nc
rm -rfv debug
