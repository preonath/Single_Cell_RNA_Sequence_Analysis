#Path declaration

sample_sheet_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Single_Cell_BCL_to_Fastq/Sample_Sheet_bcl_to_fastq"
run_folder_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Orginal_Run_Folder_NextSeq2k"
demultiplexed_path="/media/chrf/Home04/NextSeq2k/NextSeq2K_Run/Single_Cell_BCL_to_Fastq/Demultiplexed"

#Create output folder
mkdir -p -v $demultiplexed_path\/Test


#Main code
bcl2fastq --sample-sheet $sample_sheet_path\/sample_sheet_with_bcl_to_fastq_i5_only.csv -R $run_folder_path\/240109_VH00293_34_AAF5HWTM5 -o $demultiplexed_path\/Test --barcode-mismatches  1
