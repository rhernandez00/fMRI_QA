# runQualityCheck.m
 
The script reads the files found in the "input" folder and outputs files in a "txt" folder. 
The script uses movement files (.par .txt) to calculate framewise displacement or nifti files (.nii or .nii.gz) to calculate DVARS. If the script finds in the "input" folder movement and nifti files with the same name, it will calculate both measurements.

Note. If the movement file or the nifti file is not found, it will fill with zeros the corresponding column.
