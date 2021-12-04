trainingdir=$1 #'../training_data_b1000/'
newmodeldir=$2 #'testmodel'

if [ "$#" -lt 2 ]
then
	 echo "This script can be used to incrementally train UNet. If this is the first time running this script for a new model, a new config file will be generated. Otherwise training will resume from the last iteration."
	 echo ""
	 echo "Usage: $0 <directory_of_training_data> <output_directory> [optional arguments]"
	 echo ""
	 echo " -b bootstrap existing model"
	 echo " -i number of new iterations (default 100k)"
	 echo ""

	 exit 1
 fi

shift 2
iterations=100000

while getopts "b:" options; do
 case $options in
  b ) echo "bootstrapping model from $OPTARG"
	  bootstrapmodel=$OPTARG;;
  i ) echo "number of final iterations (after bootstrapping if included) $OPTARG"
	  iterations=$OPTARG;;
    * ) usage
	exit 1;;
 esac
done

if [ -f "$newmodeldir/config.ini" ]
then
python write_config_NiftyNet.py $trainingdir $newmodeldir $iterations $bootstrapmodel
else
mv $newmodeldir/dataset_split_training.csv $newmodeldir/dataset_split.csv # resume past dataset_split.csv
fi

# TO BE RUN IN BASH
# requires niftynet
net_segment -c $newmodeldir/config.ini train
net_segment -c $newmodeldir/config.ini inference
net_segment -c $newmodeldir/config.ini evaluation

# need to rename this file before AutoTops_transformAndRollOut.m
mv $newmodeldir/dataset_split.csv $newmodeldir/dataset_split_training.csv

