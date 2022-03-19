# cd /home/fs0/qijia/scratch/origin_data/cone_dev/phantom/raw_data_15-3-22
# wget -r -nd --no-parent -A ".dat" http://twix.services.fmrib.ox.ac.uk/f3t_remote/F3T_Developer/Qijia/Raw_data_15-3-22/
# wget -r -nd --no-parent -A ".pdf" http://twix.services.fmrib.ox.ac.uk/f3t_remote/F3T_Developer/Qijia/Raw_data_15-3-22/
# fsl_sub -q short.q -s openmp,8 -l logs /opt/fmrib/MATLAB/R2021a/bin/matlab -nojvm -nodisplay -batch 'recon\(1\)'
