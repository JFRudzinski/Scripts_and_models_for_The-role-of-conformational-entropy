Nproc=4
clustering="/data/isilon/rudzinski/soft_backup/Clustering-0.12_static/clustering"

R=0.1
pop_min=400
dFE=0.1

traj_len=100001
lag=20
meta_lim='0.800'


${clustering} mpp -i microstates -D fe -l ${lag} -v --concat-nframes ${traj_len}

echo "* ${lag}" &> win

${clustering} coring -s mpp_traj_${meta_lim}.dat -w win -o clustered_traj_mpp_${meta_lim} -v --concat-nframes ${traj_len}
