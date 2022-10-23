rm -rfv output.namelist*
rm -rfv run.stat*
rm -rfv ocean.output
rm -rfv DYFAMED_*_restart*.nc
toto=`awk '/output.path = /{print $NF}' output.conf`
rm -rfv $toto
rm -rfv DYFAMED_*grid_*.nc
rm -rfv DYFAMED_*ptrc_T*.nc
rm -rfv debug
