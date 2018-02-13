Nproc=16
clustering="/data/isilon/rudzinski/soft_backup/Clustering-0.12_static/clustering"

#ln -s dtraj_tica5D.dat coords 

#&& 

#${clustering} density -f coords -R 0.025 0.05 0.075 0.1 -p pop -d fe -n ${Nproc} -v 

R=0.1
pop_min=400
dFE=0.1

#&& 

ln -s fe_0.100000 fe 

#&& 

${clustering} density -f coords -r ${R} -D fe -b nn -n ${Nproc} -v 

#&& 

${clustering} density -f coords -r ${R} -D fe -B nn -T 0. $dFE -o clust -n ${Nproc} -v 

${clustering} network -p ${pop_min} --step $dFE -v && ${clustering} density -f coords -i network_end_node_traj.dat -r ${R} -D fe -B nn -o microstates -n ${Nproc} -v 

# && clustering mpp -i microstates -D fe -l 7 -v --concat-nframes 50001
