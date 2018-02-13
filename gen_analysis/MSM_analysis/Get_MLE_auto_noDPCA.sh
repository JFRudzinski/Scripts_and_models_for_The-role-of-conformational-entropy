
cd DPCA_fromTICA-5dim/

#cp ../../../../../gen_analysis/MSM_analysis/DPCA_fromTICA-5dim/commands_core_microstates.sh ./
#lag=10
#sed -i "/lag=/c\lag=${lag}" ./commands_core_microstates.sh
#./commands_core_microstates.sh

cd ../
cp -r ../../../../gen_analysis/MSM_analysis/MLE ./
cd MLE
python Convert_DPCA_dtraj_direct.py
sed -i "/lags = /c\lags = np.linspace(1,100,50).astype(int)" ./Get_MLE.py
lag=10
sed -i "/tau_CG = /c\tau_CG = ${lag}" ./Get_MLE.py
python Get_MLE.py

cd ../DPCA_fromTICA-5dim/
#cp ../../../../../gen_analysis/MSM_analysis/DPCA_fromTICA-5dim/commands_mpp.sh ./
#lag=20
#sed -i "/lag=/c\lag=${lag}" ./commands_mpp.sh
meta_lim='0.800'
#sed -i "/meta_lim=/c\meta_lim='${meta_lim}'" ./commands_mpp.sh
#./commands_mpp.sh

cd ../MLE
sed -i "/meta_lim=/c\meta_lim='${meta_lim}'" ./Convert_DPCA_dtraj_direct_mpp.py
python Convert_DPCA_dtraj_direct_mpp.py
sed -i "/lags = /c\lags = np.linspace(1,200,50).astype(int)" ./Get_MLE_mpp.py
lag=20
sed -i "/tau_CG = /c\tau_CG = ${lag}" ./Get_MLE_mpp.py
python Get_MLE_mpp.py

