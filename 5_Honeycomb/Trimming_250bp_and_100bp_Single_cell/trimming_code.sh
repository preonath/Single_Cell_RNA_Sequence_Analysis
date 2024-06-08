path="/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/Honeycomb_Biotechnologies/Data_Analysis/3_Additional_Analysis_Resources/1_Example_Data_for_BeeNetPLUS/1_Practice_running_BeeNetPLUS_using_this_dataset/FASTQs"

#trimmomatic="/home/chrf/programms/Trimmomatic-0.39/trimmomatic-0.39.jar"
#truseq3PE="/home/chrf/programms/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

for dir in $path/SC*
do 
 	name=`basename $dir  | cut -f 1,2,3,4 -d '_'`
 	echo $name
 	mkdir -p $dir\/2_Trimmed_50/
 	trimmomatic PE -threads 92 $dir\/1_RawData/$name\_*_R1_001.fastq.gz $dir\/1_RawData/$name\_*_R2_001.fastq.gz $dir\/2_Trimmed_50/$name\_R1_TrmP.fastq.gz /dev/null $dir\/2_Trimmed_50/$name\_R2_TrmP.fastq.gz /dev/null ILLUMINACLIP:$truseq3PE:2:30:10 CROP:25

done



# path="/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/Honeycomb_Biotechnologies/Data_Analysis/3_Additional_Analysis_Resources/1_Example_Data_for_BeeNetPLUS/1_Practice_running_BeeNetPLUS_using_this_dataset/FASTQs"

# #trimmomatic="/home/chrf/programms/Trimmomatic-0.39/trimmomatic-0.39.jar"
# #truseq3PE="/home/chrf/programms/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# for dir in $path/SC*
# do 
#  	name=`basename $dir  | cut -f 1,2,3,4 -d '_'`
#  	echo $name
#  	mkdir -p $dir\/2_Trimmed_50/
#  	trimmomatic SE -threads 12 $dir\/1_RawData/$name\_*_R2_001.fastq.gz  $dir\/2_Trimmed_50/$name\_R2_TrmP.fastq.gz ILLUMINACLIP:$truseq3PE:2:30:10 CROP:50

# done
