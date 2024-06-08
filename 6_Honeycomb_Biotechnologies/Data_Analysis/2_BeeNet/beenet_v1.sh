ref_path="/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/Honeycomb_Biotechnologies/Data_Analysis/3_Additional_Analysis_Resources/1_Example_Data_for_BeeNetPLUS/1_Practice_running_BeeNetPLUS_using_this_dataset"

data_path="/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/Honeycomb_Biotechnologies/Data_Analysis/3_Additional_Analysis_Resources/1_Example_Data_for_BeeNetPLUS/1_Practice_running_BeeNetPLUS_using_this_dataset/FASTQs/SC_NP_H_7061/1_RawData"

./beenet analyze --sample-name=SC_NP_H_7061 --out=SC_NP_H_7061 --num-barcodes=5000 --ref=$ref_path/20210603_GRCh38_104 $data_path/SC_NP_H_7061_S1_L001_R1_001.fastq.gz $data_path/SC_NP_H_7061_S1_L001_R2_001.fastq.gz


