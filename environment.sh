# Project-related environment variables 

#This part should be changed by the user(s)
export jobid=NOTSET
if [ -e jobid.sh ] 
then
		source jobid.sh
else
		echo "I did not find jobid.sh, are you sure you do not want to set the jobid and leave it to $jobid?"
fi
export URA_NTHREADS=1
#export URA_PROJECT_LIBS='-lwhatever_you_need_to_make_it_compile' #<-- Add here needed libraries for the project compilation

#HERE ARE LIONS!
#This part should be handled automatically by the scripts,
#touch it carefully
export URA_PROJECT=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)
#source URA env
ura_dir=`dirname $URA_PROJECT`
source $ura_dir/URAnalysis/environment.sh
export PATH=$URA_PROJECT/bin:$PATH

#FNAL PATCH
if [[ $HOSTNAME == *"fnal.gov"* ]]
then
  export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/$SCRAM_ARCH/external/boost/1.63.0/lib/:$LD_LIBRARY_PATH
  #export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/lhapdf6/6.1.5/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/$SCRAM_ARCH/external/lhapdf/6.1.6/lib:$LD_LIBRARY_PATH
fi

export ENV_SOURCED="yes"

export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH

#export LHAPDF_DATA_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/lhapdf6/6.1.5/share/LHAPDF
