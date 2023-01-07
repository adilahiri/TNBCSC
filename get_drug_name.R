get_drug_name <- function(mat,num_drug,n){

  get_index <- which(mat[,13]==num_drug, arr.ind = TRUE)
  drug_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(drug_df)<- c("Drugs Combination", "Size Difference")
  for (iter in 1:n){
   
    if(num_drug>1){
    col_index <- which(mat[get_index[iter],]==1)  
    drug_df[iter,1]=paste(colnames(mat)[col_index],collapse=" + ")
    drug_df[iter,2]=mat[get_index[iter],12]
    }
    else
    {
      col_index <- which(mat[get_index[iter],]==1)[1]
      drug_df[iter,1]=colnames(mat)[col_index]
      drug_df[iter,2]=mat[get_index[iter],12]
    }
  }
 return(drug_df)
}