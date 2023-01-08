source('~/Desktop/Breast_Cancer/code/TNBC_one_fault.R')
source('~/Desktop/Breast_Cancer/code/TNBC_two_fault.R')
source('~/Desktop/Breast_Cancer/code/TNBC_three_fault.R')
source('~/Desktop/Breast_Cancer/code/get_drug_name.R')


library(openxlsx)
library(dplyr)
library(xlsx)
library(ggplot2)
#require(gridExtra)
library(patchwork)

print("Starting anlysis")

start_time <- Sys.time() # get start time
# print elapsed time
new <- Sys.time() - start_time 
print("time elapsed:") 
print(new)

 
N= 5 # number of drugs
A= matrix(0,43,2048) # store SD matrix
Drug= matrix(0,1,11) # map for drug vector
D_iter=1

drug_name <- c("BCAT_comple_Inhibitor","GLI_Inhibitor","Notch_complex_Inhibitor",
               "AKT_Inhibitor","NFKB_Inhibitor","NIK_Inhibitor","CREB_Inhibitor",
               "STAT3_Inhibitor","PI3K_Inhibitor","SMAD_complex_Inhibitor",
               "Hippo_complex_Inhibitor")
print("Computing SD matrix for one fault at a time")
for (i in 0:42){
  for (d1 in 0:1)
    for(d2 in 0:1)
      for(d3 in 0:1)
        for(d4 in 0:1)
          for(d5 in 0:1)
            for(d6 in 0:1)
              for(d7 in 0:1)
                for(d8 in 0:1)
                  for(d9 in 0:1)
                    for(d10 in 0:1)
                      for(d11 in 0:1){
                        m = 1024*d1 + 512*d2 + 256*d3 + 128*d4 + 64*d5 + 32*d6 + 16*d7 + 8*d8 + 4*d9+ 2*d10+ d11+1
                        if (d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11 <= N) {
                          A[i+1,m] = TNBC_one_fault(i,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11);
                          drug_vector=c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11);
                          Drug= rbind(Drug,drug_vector);
                          D_iter=D_iter+1;
                        }
                      }
}

Drug <- Drug [-1,]
A1=colSums(A)
indices<-which(A1==0)
A1<-A1[-indices]
A1<-A1/max(abs(A1))
one_mutation_df<- as.data.frame(cbind(Drug[1:1024,],A1))
colnames(one_mutation_df)<-c(drug_name,"size_difference")
one_mutation_df<- one_mutation_df %>% mutate(number_of_drugs_applied = rowSums(.[1:11]))
print ("Completed computing SD matrix for one fault at a time")
# print elapsed time
new <- Sys.time() - start_time 
print("time elapsed:") 
print(new)

# Two Fault
B<-array(0, dim = c(43, 43, 2048))
print("Computing SD matrix for two fault a time")
for (i in 0:42){
 for(j in 0:42)  
    for (d1 in 0:1)
      for(d2 in 0:1)
        for(d3 in 0:1)
          for(d4 in 0:1)
            for(d5 in 0:1)
              for(d6 in 0:1)
                for(d7 in 0:1)
                  for(d8 in 0:1)
                    for(d9 in 0:1)
                      for(d10 in 0:1)
                        for(d11 in 0:1){
                          m = 1024*d1 + 512*d2 + 256*d3 + 128*d4 + 64*d5 + 32*d6 + 16*d7 + 8*d8 + 4*d9+ 2*d10+ d11+1
                          if (d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11 <= N) {
                            B[i+1,j+1,m] = TNBC_two_fault(i,j,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11);
                     
                          }
                      }
}


B1=colSums(colSums(B))
indices<-which(B1==0)
B1<-B1[-indices]
B1<-B1/max(abs(B1))

two_mutation_df<- as.data.frame(cbind(Drug[1:1024,],B1))
colnames(two_mutation_df)<-c(drug_name,"size_difference")
two_mutation_df<- two_mutation_df %>% mutate(number_of_drugs_applied = rowSums(.[1:11]))
print ("Completed computing SD matrix for two fault at a time")
# print elapsed time
new <- Sys.time() - start_time 
print("time elapsed:") 
print(new)
#Three Fault Network
C<-array(0, dim = c(43,43, 43, 2048))

print("Computing SD matrix for three fault a time")

for (i in 0:42){
  for(j in 0:42)
    for(k in 0:42)
      for (d1 in 0:1)
        for(d2 in 0:1)
          for(d3 in 0:1)
            for(d4 in 0:1)
              for(d5 in 0:1)
                for(d6 in 0:1)
                  for(d7 in 0:1)
                    for(d8 in 0:1)
                      for(d9 in 0:1)
                        for(d10 in 0:1)
                          for(d11 in 0:1){
                            m = 1024*d1 + 512*d2 + 256*d3 + 128*d4 + 64*d5 + 32*d6 + 16*d7 + 8*d8 + 4*d9+ 2*d10+ d11+1
                            if (d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11 <= N) {
                              C[i+1,j+1,k+1,m] = TNBC_three_fault(i,j,k,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11);
                          }
                    }
}

C1=colSums(colSums(colSums(C)))
indices<-which(C1==0)
C1<-C1[-indices]
C1<-C1/max(abs(C1))

three_mutation_df<- as.data.frame(cbind(Drug[1:1024,],C1))
colnames(three_mutation_df)<-c(drug_name,"size_difference")
three_mutation_df<- three_mutation_df %>% mutate(number_of_drugs_applied = rowSums(.[1:11]))
print ("Completed computing SD matrix for three fault at a time")
# print elapsed time
new <- Sys.time() - start_time 
print("time elapsed") 
print(new)

# Generate tables
one_mutation_df <- one_mutation_df %>% arrange(size_difference)
two_mutation_df <- two_mutation_df %>% arrange(size_difference)
three_mutation_df <- three_mutation_df %>% arrange(size_difference)



# Generate top 10 scores for single mutation network with different combination of drugs
top10_one_drug_one_mutation <- get_drug_name(one_mutation_df,num_drug = 1,n=10)
top10_two_drug_one_mutation <- get_drug_name(one_mutation_df,num_drug = 2,n=10)
top10_three_drug_one_mutation <- get_drug_name(one_mutation_df,num_drug = 3,n=10)
top10_four_drug_one_mutation <- get_drug_name(one_mutation_df,num_drug = 4,n=10)
top10_five_drug_one_mutation <- get_drug_name(one_mutation_df,num_drug = 5,n=10)


one_mutation_names <- list('All_Drug_Combos' = one_mutation_df, 'One_Drug' = top10_one_drug_one_mutation, 
                           'Two_Drug' = top10_two_drug_one_mutation,'Three_Drug' = top10_three_drug_one_mutation,
                           'Four_Drug' = top10_four_drug_one_mutation,'Five_Drug' = top10_five_drug_one_mutation)
write.xlsx(one_mutation_names, file = 'One_Mutation_File.xlsx')


# Generate top 10 scores for double mutation network with different combination of drugs
top10_one_drug_two_mutation <- get_drug_name(two_mutation_df,num_drug = 1,n=10)
top10_two_drug_two_mutation <- get_drug_name(two_mutation_df,num_drug = 2,n=10)
top10_three_drug_two_mutation <- get_drug_name(two_mutation_df,num_drug = 3,n=10)
top10_four_drug_two_mutation <- get_drug_name(two_mutation_df,num_drug = 4,n=10)
top10_five_drug_two_mutation <- get_drug_name(two_mutation_df,num_drug = 5,n=10)


two_mutation_names <- list('All_Drug_Combos' = two_mutation_df, 'One_Drug' = top10_one_drug_two_mutation, 
                           'Two_Drug' = top10_two_drug_two_mutation,'Three_Drug' = top10_three_drug_two_mutation,
                           'Four_Drug' = top10_four_drug_two_mutation,'Five_Drug' = top10_five_drug_two_mutation)
write.xlsx(two_mutation_names, file = 'Two_Mutation_File.xlsx')

# Generate top 10 scores for triple mutation network with different combination of drugs
top10_one_drug_three_mutation <- get_drug_name(three_mutation_df,num_drug = 1,n=10)
top10_two_drug_three_mutation <- get_drug_name(three_mutation_df,num_drug = 2,n=10)
top10_three_drug_three_mutation <- get_drug_name(three_mutation_df,num_drug = 3,n=10)
top10_four_drug_three_mutation <- get_drug_name(three_mutation_df,num_drug = 4,n=10)
top10_five_drug_three_mutation <- get_drug_name(three_mutation_df,num_drug = 5,n=10)


three_mutation_names <- list('All_Drug_Combos' = three_mutation_df, 'One_Drug' = top10_one_drug_three_mutation, 
                           'Two_Drug' = top10_two_drug_three_mutation,'Three_Drug' = top10_three_drug_three_mutation,
                           'Four_Drug' = top10_four_drug_three_mutation,'Five_Drug' = top10_five_drug_three_mutation)
write.xlsx(three_mutation_names, file = 'Three_Mutation_File.xlsx')

## Make plots 
top10_one_drug_one_mutation<- rbind(top10_one_drug_one_mutation,c("Untreated",1))
top10_two_drug_one_mutation<- rbind(top10_two_drug_one_mutation,c("Untreated",1))
top10_three_drug_one_mutation<- rbind(top10_three_drug_one_mutation,c("Untreated",1))
top10_four_drug_one_mutation<- rbind(top10_four_drug_one_mutation,c("Untreated",1))
top10_five_drug_one_mutation<- rbind(top10_five_drug_one_mutation,c("Untreated",1))

top10_one_drug_two_mutation <- rbind(top10_one_drug_two_mutation,c("Untreated",1))
top10_two_drug_two_mutation <- rbind(top10_two_drug_two_mutation,c("Untreated",1))
top10_three_drug_two_mutation <-rbind(top10_three_drug_two_mutation,c("Untreated",1))
top10_four_drug_two_mutation <- rbind(top10_four_drug_two_mutation,c("Untreated",1))
top10_five_drug_two_mutation <- rbind(top10_five_drug_two_mutation,c("Untreated",1))

top10_one_drug_three_mutation <- rbind(top10_one_drug_three_mutation,c("Untreated",1))
top10_two_drug_three_mutation <- rbind(top10_two_drug_three_mutation,c("Untreated",1))
top10_three_drug_three_mutation <- rbind(top10_three_drug_three_mutation,c("Untreated",1))
top10_four_drug_three_mutation <- rbind(top10_four_drug_three_mutation,c("Untreated",1))
top10_five_drug_three_mutation <- rbind(top10_five_drug_three_mutation,c("Untreated",1))



# Plot of single mutations 
one_mut_one_drug_pt<-ggplot(top10_one_drug_one_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

one_mut_two_drug_pt<-ggplot(top10_two_drug_one_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

one_mut_three_drug_pt<-ggplot(top10_three_drug_one_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

one_mut_four_drug_pt<-ggplot(top10_four_drug_one_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

one_mut_five_drug_pt<-ggplot(top10_five_drug_one_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 


P1<- one_mut_one_drug_pt + one_mut_two_drug_pt + one_mut_three_drug_pt +
  one_mut_four_drug_pt + one_mut_five_drug_pt + plot_annotation(tag_levels = "I") 

#grid.arrange(one_mut_one_drug_pt, one_mut_two_drug_pt, one_mut_three_drug_pt,
#             one_mut_four_drug_pt,one_mut_five_drug_pt, ncol=3,nrow=2)


# Plot of double mutations 

two_mut_one_drug_pt<-ggplot(top10_one_drug_two_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

two_mut_two_drug_pt<-ggplot(top10_two_drug_two_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

two_mut_three_drug_pt<-ggplot(top10_three_drug_two_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

two_mut_four_drug_pt<-ggplot(top10_four_drug_two_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

two_mut_five_drug_pt<-ggplot(top10_five_drug_two_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 


P2<- two_mut_one_drug_pt + two_mut_two_drug_pt + two_mut_three_drug_pt + 
  two_mut_four_drug_pt + two_mut_five_drug_pt + plot_annotation(tag_levels = "I")


# Plot of triple mutations 

three_mut_one_drug_pt<-ggplot(top10_one_drug_three_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

three_mut_two_drug_pt<-ggplot(top10_two_drug_three_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

three_mut_three_drug_pt<-ggplot(top10_three_drug_three_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

three_mut_four_drug_pt<-ggplot(top10_four_drug_three_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 

three_mut_five_drug_pt<-ggplot(top10_five_drug_three_mutation, aes(x=Drugs_Combination, y=Norm_Mean_Size_Difference)) +geom_bar(stat="identity")+
  theme(text = element_text(size=5),axis.text.x = element_text(angle=45, hjust=1)) 


P3<-three_mut_one_drug_pt + three_mut_two_drug_pt + three_mut_three_drug_pt + 
  three_mut_four_drug_pt +three_mut_five_drug_pt + plot_annotation(tag_levels = "I")


ggsave("Single_Mutation_Drug_SD.png",plot=P1,dpi=300)
ggsave("Double_Mutation_Drug_SD.png",plot=P1,dpi=300)
ggsave("Triple_Mutation_Drug_SD.png",plot=P1,dpi=300)

print("Analysis complete, tables generated in same directory")


