Nproc=4
clustering="/data/isilon/rudzinski/soft_backup/Clustering-0.12_static/clustering"

R=0.1
pop_min=400
dFE=0.1

traj_len=100001
lag=10
meta_lim='0.750'


echo "* ${lag}" &> win

${clustering} coring -s microstates -w win -o microstates_cored -v --concat-nframes ${traj_len}
