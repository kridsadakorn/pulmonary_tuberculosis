ipcaps cleaned_data.bed -o pca_results

#For loop files, you can use shell script :
for((k=2;k<21;k++)); do admixture cleaned_data.bed k --cv=10 -j4; done
#Where :
# k = The number of populations (K) that was assumed for the analysis.
# --cv = Cross-validation, 10-fold CV was used in this case. The cross-validation error is reported in the output.
# -j = The number of thread
