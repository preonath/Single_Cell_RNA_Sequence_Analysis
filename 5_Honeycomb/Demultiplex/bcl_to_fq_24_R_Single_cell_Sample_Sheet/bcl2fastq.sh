path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/bcl_to_fq_24_R_Single_cell_Sample_Sheet"     
bcl2fastq --sample-sheet $path/index_1_Nor_and_index_2_F_Nor.csv -R 240109_VH00293_34_AAF5HWTM5 -o index_1_Nor_and_index_2_F_Nor_Fastq_NS2k_Batch_24_R_SC/ --barcode-mismatches 2

path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/bcl_to_fq_24_R_Single_cell_Sample_Sheet" 
bcl2fastq --sample-sheet $path/index_1_Nor_and_index_2_F_Rev_com.csv -R 240109_VH00293_34_AAF5HWTM5 -o index_1_Nor_and_index_2_F_Rev_com_Fastq_NS2k_Batch_24_R_SC_mismatch_2/ --barcode-mismatches  2


path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/bcl_to_fq_24_R_Single_cell_Sample_Sheet" 
bcl2fastq --sample-sheet $path/index_1_Rev_com_and_index_2_F_Nor.csv -R 240109_VH00293_34_AAF5HWTM5 -o index_1_Rev_com_and_index_2_F_Nor_Fastq_NS2k_Batch_24_R_SC/ --barcode-mismatches  2


path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/bcl_to_fq_24_R_Single_cell_Sample_Sheet" 
bcl2fastq --sample-sheet $path/index_1_Rev_com_and_index_2_F_Rev_com.csv -R 240109_VH00293_34_AAF5HWTM5 -o index_1_Rev_com_and_index_2_F_Rev_com_Fastq_NS2k_Batch_24_R_SC/ --barcode-mismatches  1

