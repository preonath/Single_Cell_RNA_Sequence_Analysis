data_path_1="/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/Single_Cell_NP_NextSeq_Batch_24_R/SC_NP_H_7061/1_RawData"
data_path_2="/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/Single_Cell_NP_NextSeq_Batch_24_R/SC_NP_H_7063/1_RawData"
ref_path="/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/BeeNet/20210603_GRCh38.104"

beenet analyze --sample-name=SC_NP_H_7061  --num-barcodes=6000 --ref=$ref_path $data_path_1/SC_NP_H_7061_S1_L001_R1_001.fastq.gz $data_path_1/SC_NP_H_7061_S1_L001_R2_001.fastq.gz $data_path_2/SC_NP_H_7063_S2_L001_R1_001.fastq.gz $data_path_2/SC_NP_H_7063_S2_L001_R1_001.fastq.gz



data_path_1="/media/chrf/Home04/NextSeq2k/Preonath_Data_Share/Single_cell/Single_Cell_NP_NextSeq_Batch_24_R/SC_NP_H_7061/1_RawData"
data_path_2="/media/chrf/Home04/NextSeq2k/Preonath_Data_Share/Single_cell/Single_Cell_NP_NextSeq_Batch_24_R/SC_NP_H_7063/1_RawData"
ref_path="/media/chrf/Home03/Trial_Preonath/Single_Cell/20210603_GRCh38.104"

beenet analyze --sample-name=SC_NP_H_7061  --num-barcodes=6000 --ref=$ref_path $data_path_1/SC_NP_H_7061_S1_L001_R1_001.fastq.gz $data_path_1/SC_NP_H_7061_S1_L001_R2_001.fastq.gz $data_path_2/SC_NP_H_7063_S2_L001_R1_001.fastq.gz $data_path_2/SC_NP_H_7063_S2_L001_R1_001.fastq.gz

/media/chrf/Home04/NextSeq2k/Preonath_Data_Share/Single_cell/Single_Cell_NP_NextSeq_Batch_24_R
