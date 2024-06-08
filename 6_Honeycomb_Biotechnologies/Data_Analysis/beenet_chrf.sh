beenet_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Honeycomb_Biotechnologies/Data_Analysis"

ref_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Honeycomb_Biotechnologies/Data_Analysis"

data_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Single_Cell_BCL_to_Fastq/Demultiplexed/index_1_Nor_and_index_2_F_Rev_com_Fastq_NS2k_Batch_24_R_SC_i5_only_v2"

result_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Honeycomb_Biotechnologies/Data_Analysis/Analysis_Result/Beenet_SC_NP_H_100"

$beenet_path/beenet analyze --sample-name=Beenet_SC_NP_H_100 --out=$result_path/Beenet_SC_NP_H_100 --num-barcodes=6000 --ref=$ref_path/20210603_GRCh38_104 $data_path/SC_NP_H_7061/2_Trimmed/SC_NP_H_7061_R1_TrmP.fastq.gz $data_path/SC_NP_H_7061/2_Trimmed/SC_NP_H_7061_R2_TrmP.fastq.gz $data_path/SC_NP_H_7063/2_Trimmed/SC_NP_H_7063_R1_TrmP.fastq.gz $data_path/SC_NP_H_7063/2_Trimmed/SC_NP_H_7063_R1_TrmP.fastq.gz


