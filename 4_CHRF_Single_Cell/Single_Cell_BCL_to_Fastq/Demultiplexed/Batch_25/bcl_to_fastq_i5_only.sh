#Path declaration

sample_sheet_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Single_Cell/Single_Cell_BCL_to_Fastq/Demultiplexed/Batch_25"
run_folder_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Orginal_Run_Folder_NextSeq2k"
demultiplexed_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Single_Cell/Single_Cell_BCL_to_Fastq/Demultiplexed/Batch_25/Demultiplex"

#Create output folder
mkdir -p -v $demultiplexed_path\/Test


#Main code
bcl2fastq --sample-sheet $sample_sheet_path\/NextSeq_Batch_25_P2_Sample_Sheet_local_29Feb2024_bcl_to_fq.csv -R $run_folder_path\/240229_VH00293_35_AACVNMKM5 -o $demultiplexed_path/Test --barcode-mismatches  1
