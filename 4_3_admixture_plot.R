
##########################
# Create a plot for cross-validation error of ADMIXTURE's results

# Read the file
file_path = "result_ADMIXTURE/result_CV_error.txt"
cv_data = readLines(file_path)

# Initialize empty vectors to store K values and CV errors
K_values = numeric()
CV_errors = numeric()

# Loop through each line and extract K and CV error values
for (line in cv_data) {
  # Extract K value
  K_value = as.numeric(sub(".*\\(K=(\\d+)\\).*", "\\1", line))
  
  # Extract CV error value
  CV_error = as.numeric(sub(".*: (\\d+\\.\\d+)$", "\\1", line))
  
  # Append values to the vectors
  K_values = c(K_values, K_value)
  CV_errors = c(CV_errors, CV_error)
}

# Sort the vectors by K values
sorted_indices = order(K_values)
K_values_sorted = K_values[sorted_indices]
CV_errors_sorted = CV_errors[sorted_indices]

# Create the desired vectors
xlabel = K_values_sorted
CV = CV_errors_sorted

xlable.at = 1:length(CV)
xlabel.txt = as.character(xlabel)

# Create a plot and save as PDF
fname="result_ADMIXTURE/admixture_k_cv.pdf"
pdf(fname)
plot(CV,type="l",xlab = "K",ylab="CV" ,xaxt="n")
axis(side=1,at = xlable.at, labels = xlabel.txt)
dev.off()

##########################
# Create a plot for ADMIXTURE's results 

ipcaps.result="result_IPCAPS"
admixture.result="result_ADMIXTURE/model.jody.14."

groupfile=paste0(ipcaps.result,"/groups.txt")
groups = read.table(file=groupfile,header=T)

# Create the table 
table_result = table(filter_group[, c('node', 'label')])

# Convert the table to a dataframe 
df_result = as.data.frame(table_result) 
# Save the dataframe to a CSV file 
write.csv(df_result, "result_IPCAPS/table_result.csv", row.names = FALSE)


node = table(groups$node)
idx = which(node >= 10)

# Manually sort the plots, make it beautiful
node.order=c(35, 60, 70, 71, 23, 27, 32, 29, 72, 38, 54, 67, 66, 20)

myspace = c()
for (i in node.order){
  ln = length(groups$row.number[which(groups$node == i)])
  tmp = rep(0,(ln-1))
  myspace = c(myspace,5,tmp)
}

# Create a plot for K = 3

fname = paste0(admixture.result,"3.Q")
order.group = NULL
my.mat = read.table(fname)
my.mat = as.matrix(my.mat)[,c(2, 1, 3)]
for (i in node.order){
  tmp = groups$row.number[which(groups$node == i)]
  print(i)
  tmp.mat = my.mat[tmp,]
  if (length(tmp)>1){
    tmp.mat.order = NULL
    tmp.mat.order = tmp.mat[order(tmp.mat[,1],tmp.mat[,2],tmp.mat[,3],decreasing=FALSE),]
    order.group = rbind(order.group,tmp.mat.order)
  }else{
    order.group = rbind(order.group,tmp.mat)
  }
}
order.group.k3 = t(order.group)

# Create a plot for K = 3

fname = paste0(admixture.result,"4.Q")
order.group = NULL
my.mat = read.table(fname)
my.mat = as.matrix(my.mat)[,c(3,1,4,2)]
for (i in node.order){
  tmp = groups$row.number[which(groups$node == i)]
  tmp.mat = my.mat[tmp,]
  if (length(tmp)>1){
    tmp.mat.order = NULL
    tmp.mat.order = tmp.mat[order(tmp.mat[,1],tmp.mat[,2],tmp.mat[,3],tmp.mat[,4],decreasing=FALSE),]
    order.group = rbind(order.group,tmp.mat.order)
  }else{
    order.group = rbind(order.group,tmp.mat)
  }
}
order.group.k4 = t(order.group)

# Create a plot for K = 5

fname = paste0(admixture.result,"5.Q")
order.group = NULL
my.mat = read.table(fname)
my.mat = as.matrix(my.mat)[,c(2,5,4,1,3)]
for (i in node.order){
  tmp = groups$row.number[which(groups$node == i)]
  tmp.mat = my.mat[tmp,]
  if (length(tmp)>1){
    tmp.mat.order = NULL
    tmp.mat.order = tmp.mat[order(tmp.mat[,1],tmp.mat[,2],tmp.mat[,3],tmp.mat[,4],tmp.mat[,5],decreasing=FALSE),]
    order.group = rbind(order.group,tmp.mat.order)
  }else{
    order.group = rbind(order.group,tmp.mat)
  }
}
order.group.k5 = t(order.group)

# Save a plot to PDF

fname="result_ADMIXTURE/admixture_profile.pdf"
pdf(fname)
par(mfrow=c(3,1))

c3 = c("#00FFFFFF","#8000FFFF","#FDE725FF")
c4 = c("#7AD151FF","#00FFFFFF","#8000FFFF","#FDE725FF")
c5 = c("#FF0000FF","#7AD151FF","#00FFFFFF","#8000FFFF","#FDE725FF")

barplot(order.group.k3, col=c3, xlab="", ylab="Ancestry (K=3)", border=NA, space=myspace)
barplot(order.group.k4, col=c4, xlab="", ylab="Ancestry (K=4)", border=NA, space=myspace)
barplot(order.group.k5, col=c5, xlab="", ylab="Ancestry (K=5)", border=NA, space=myspace)
dev.off()

#################
# Create a plot for PCA 

library(readxl)
library(KRIS)

ipcaps.result="selected_model/IPCAPs_result"
#groupfile=paste0(ipcaps.result,"/groups_3nodes.xls")
#groups = read_excel(groupfile, sheet = "groups_3nodes")
groupfile=paste0(ipcaps.result,"/groups.txt")
groups = read.delim(groupfile)

nodefile = paste0(ipcaps.result,"/RData/node1.RData")
load(nodefile)

node.order=c(35, 60, 70, 71, 23, 27, 32, 29, 72, 38, 54, 67, 66, 20)
new_group = groups[order(groups$row.number),]

idx = which(new_group$node %in% node.order)
filter_group = new_group[idx,]
plot_PC = PCs[filter_group$row.number,]
plot_label = filter_group$node

fname="selected_model/PC_14nodes.pdf"
pdf(fname, width = 10, height = 10)

my.pattern = rep(1,14)


plot_label[which(plot_label == 35)] = 'Node 35'
plot_label[which(plot_label == 60)] = 'Node 60'
plot_label[which(plot_label == 70)] = 'Node 70'
plot_label[which(plot_label == 71)] = 'Node 71'
plot_label[which(plot_label == 23)] = 'Node 23'
plot_label[which(plot_label == 27)] = 'Node 27'
plot_label[which(plot_label == 32)] = 'Node 32'
plot_label[which(plot_label == 29)] = 'Node 29'
plot_label[which(plot_label == 72)] = 'Node 72'
plot_label[which(plot_label == 38)] = 'Node 38'
plot_label[which(plot_label == 54)] = 'Node 54'
plot_label[which(plot_label == 67)] = 'Node 67'
plot_label[which(plot_label == 66)] = 'Node 66'
plot_label[which(plot_label == 20)] = 'Node 20'
all.labels = unique(plot_label)

map.color = c("#0000FF","#00FF00","#FF0000","#00FFFF","#FF00FF",
              "#FFFF00","#008888","#880088","#888800","#000088",
              "darkorange","#880000","#888888", "#000000")

plot3views(plot_PC, labels = plot_label, plot.legend = all.labels, plot.pattern = my.pattern, plot.color = map.color)

dev.off()


###############################
# Calculate Fst values for clustering results

library(KRIS)

ipcaps.result="selected_model/IPCAPs_result"
# load raw data
rawfile = paste0(ipcaps.result,"/RData/rawdata.RData")
load(rawfile)

groupfile=paste0(ipcaps.result,"/groups.txt")
groups = read.delim(groupfile)

tmp = table(groups$node)
valid_node = as.factor(names(tmp[tmp>10]))

n = length(valid_node)
mat_fst = matrix(rep(0, n*n), nrow = n)
colnames(mat_fst) = valid_node
rownames(mat_fst) = valid_node

for (i in 1:(length(valid_node)-1)){
  for (j in (i+1):length(valid_node)){
    print(paste0("i = ",i,", j = ",j))
    idx1 = groups$row.number[which(groups$node == valid_node[i])]
    idx2 = groups$row.number[which(groups$node == valid_node[j])]
    v_fst = fst.hudson(raw.data, idx1, idx2)
    print(paste0("Fst = ", v_fst))
    mat_fst[i,j] = v_fst
  }
}

fout = "selected_model/pairwise_fst_all_nodes.csv"
write.csv(mat_fst, file = fout)
