TNBC_two_fault <- function(fault1, fault2, x1, x2, x3, x4, x5,x6,x7,x8,x9,x10,x11) {
  
  
  # Input Vectors in ideal state
  WNT=0; HH=0;Jag=0;Delta=0;Hippo=1;EGF=0;
  IGF=0;PTEN=1;TGFB=0;SMAD4=0;SMAD7=1;
  TNF=0;IL6=0;IL8=0;LEP=0;HNL1=0;
  CXCL1=0;CXCL3=0;HAS1=0;PTGIS=0;PFKFB3=0;
  
  # Drug_Vector
  BCAT_Inhi = x1;
  GLI_Inhi = x2;
  Notch_complex_Inhi =x3;
  AKT_Inhi =x4;
  NFKB_Inhi=x5;
  NIK_Inhi=x6;
  CREB_Inhi= x7;
  Crypto =x8;
  PI3K_Inhi= x9;
  SMAD_Inhi=x10;
  YAP_TAZ_TEAD_Inhi=x11;
  
  if (fault1 == 1 || fault2==1){
    FRZ= 1; 
  } 
  else{ 
    FRZ=WNT;
  }  
  
  
  if (fault1 == 2 || fault2==2){
    DVL= 0;    
  }
  else{
    DVL= !FRZ;
  }
  if (fault1 == 3 || fault2==3){
    GSK_Axin =1; 
  }
  else{
    GSK_Axin = !DVL;
  }
  
  if (fault1 == 4 || fault2==4){
    Beta_Cat =1;
  }
  else{
    Beta_Cat = GSK_Axin;
  }
  
  if (fault1 == 5 || fault2==5){
    TCF_LEF=1;
  }
  else {
    TCF_LEF = Beta_Cat &(!BCAT_Inhi);
  }
  if (fault1 == 6 || fault2==6) {
    HH_Inv=0;     
  }
  else{
    HH_Inv = !HH;
  }
  
  if (fault1 == 7 || fault2==7) {
    PTCH=1;    
  }
  else{
    PTCH= ! HH_Inv;
  }
  
  if (fault1 == 8 || fault2==8){
    SMO=1;     
  }
  else{
    SMO= PTCH;
  }
  
  if (fault1 == 9 || fault2==9){
    GLI=1;    
  }
  else{
    GLI = SMO;
  }
  
  if (fault1 == 10 || fault2==10)  {
    NOTCH=1;    
  }
  else{
    NOTCH = Delta & Jag;
  }
  
  if (fault1 == 11 || fault2==11){
    NICD =1;   
  }
  else{
    NICD= NOTCH;
  }
  
  if (fault1 == 12 || fault2==12) {
    CSL_MAML_RBJK=1;
  }
  else{
    CSL_MAML_RBJK = NICD;
  }
  
  if (fault1 == 13 || fault2==13){
    SAV1_MST =0;   
  }
  else{
    SAV1_MST= Hippo;
  }
  
  if (fault1 == 14 || fault2==14){
    LST_MOB =1;   
  }
  else{
    LST_MOB= ! SAV1_MST;
  }
  
  if (fault1 == 15 || fault2==15){
    YAP_TAZ=1;    
  }
  else {
    YAP_TAZ = LST_MOB;
  }  
  
  
  if (fault1 == 16 || fault2==16){
    YAP_TAZ_TEAD=1;
  }
  else{
    YAP_TAZ_TEAD =YAP_TAZ;
  }
  
  
  if (fault1 == 17 || fault2==17){
    EGFR=1;    
  }
  else{
    EGFR = EGF;
  }
  
  
  if (fault1 == 18 || fault2==18){
    IGF1R=1;    
  }
  else {
    IGF1R = IGF;
  }
  
  
  if (fault1 == 19 || fault2==19) {
    PTEN_INV=1;    
  }
  else{
    PTEN_INV= !PTEN;
  }
  
  
  if (fault1 == 20 || fault2==20) {
    TGFBR2=1;    
  }
  else{
    TGFBR2= TGFB;
  }
  
  if (fault1 == 21 || fault2==21) {
    TGFBR1=1;    
  }
  else{
    TGFBR1=TGFBR2;      
  }
  
  if (fault1 == 22 || fault2==22) {
    MAP3K7=1;    
  }
  else{
    MAP3K7= TGFBR1;
  }
  
  if (fault1 == 23 || fault2==23){
    SMAD7_INV=1;  
  }
  else{
    SMAD7_INV = !SMAD7;
  }
  
  if (fault1 == 24 || fault2==24) {
    SMAD2_3=1;    
  }
  else {
    SMAD2_3= TGFBR1 & SMAD4 & SMAD7_INV;
  }
  
  
  if (fault1 == 25 || fault2==25) {
    SMAD2_3_4=1;   
  }
  else{
    SMAD2_3_4 = SMAD2_3 & (!SMAD_Inhi);   
  }
  
  if (fault1 == 26 || fault2==26) {
    TNFR2=1;    
  }
  else {
    TNFR2= TNF;
  }
  
  
  if (fault1 == 27 || fault2==27) {
    TRAF2=1;    
  }
  else {
    TRAF2= TNFR2;
  }
  
  
  
  if (fault1 == 28 || fault2==28) {
    PI3K=1;    
  }
  else {
    PI3K = EGFR | IGF1R |TGFBR1 | TRAF2 ;
  }
  
  if (fault1 == 29 || fault2==29) {
    PIP3=1;   
  }
  else{
    PIP3= (PI3K &(!PI3K_Inhi)) | PTEN_INV;  
  }
  
  if (fault1 == 30 || fault2==30) {
    AKT=1;    
  }
  else{
    AKT= PIP3;
  }
  
  
  if (fault1 == 31 || fault2==31) {
    NIK=1;     
  }
  else{
    NIK= TRAF2;
  }
  
  
  if (fault1 == 32 || fault2==32) {
    CREB=1;    
  }
  else{
    CREB= AKT &(!AKT_Inhi);
  }
  
  
  if (fault1 == 33 || fault2==33) {
    AKT_INV=0;    
  }
  else{
    AKT_INV = !(AKT& (!AKT_Inhi));
  }
  
  
  if (fault1 == 34 || fault2==34) {
    IKK=1;     
  }
  else {
    IKK= (AKT& (!AKT_Inhi)) | (NIK &(!NIK_Inhi)); 
  }
  
  
  if (fault1 == 35 || fault2==35) {
    IKBa=1;     
  }
  else{
    IKBa= IKK;
  }
  
  
  if (fault1 == 36 || fault2==36) {
    P27=1;  
  }
  else {
    P27= !AKT_INV; 
  }
  
  if (fault1 == 37 || fault2==37) {
    IL6R=1;    
  }
  else {
    IL6R= IL6;
  }
  
  if (fault1 == 38 || fault2==38) {
    IL8R=1;    
  }
  else{
    IL8R= IL8;
  }
  
  
  if (fault1 == 39 || fault2==39) {
    LEPR=1;    
  }
  else {
    LEPR= LEP | HNL1;
  }
  
  
  if (fault1 == 40 || fault2==40) {
    JAK2=1;    
  }
  else{
    JAK2= IL6R | IL8R | LEPR |CXCL1| CXCL3|HAS1|PTGIS|PFKFB3;
  }
  
  
  if (fault1 == 41 || fault2==41) {
    STAT3=1;    
  }
  else {
    STAT3=JAK2;       
  }
  
  if (fault1 == 42 || fault2==42) {
    NFKB=1;    
  }
  else {
    NFKB= (CSL_MAML_RBJK & (!Notch_complex_Inhi)) |IKBa | MAP3K7 | (STAT3 & (!Crypto));       
  }
  
  CCND1 = TCF_LEF | (GLI & (!GLI_Inhi) ) | (CSL_MAML_RBJK & (!Notch_complex_Inhi)) | (YAP_TAZ_TEAD & (!YAP_TAZ_TEAD_Inhi)) | P27| (STAT3 & (!Crypto));
  
  C_Myc= TCF_LEF| (GLI & (!GLI_Inhi)) | (CSL_MAML_RBJK & (!Notch_complex_Inhi)) | (YAP_TAZ_TEAD &(!YAP_TAZ_TEAD_Inhi)) |(STAT3 & (!Crypto));
  
  SNAIL= TCF_LEF|(GLI & (!GLI_Inhi))| (NFKB & (!NFKB_Inhi)) |SMAD2_3_4;
  
  BCL2= (GLI & (!GLI_Inhi))|(CREB &(!CREB_Inhi))|(STAT3 & (!Crypto));
  
  TWIST=TCF_LEF|(GLI & (!GLI_Inhi))| (NFKB & (!NFKB_Inhi));
  
  MCL1=(CREB &(!CREB_Inhi))|(STAT3 & (!Crypto));
  
  gate_output <- as.numeric(c(CCND1, C_Myc, SNAIL,BCL2, TWIST,MCL1));
  ideal_output <- c(0,0,0,0,0,0);
  
  a=0; b=0; c=0; d=0;
  
  for (i in 1:6){
    if (gate_output[i] == 1 && ideal_output[i] == 1){
      a=a+1; }
    else if (gate_output[i] == 0 && ideal_output[i] == 1){
      b=b+1;} 
    else if (gate_output[i] == 1 && ideal_output[i] == 0){
      c=c+1; }
    else if (gate_output[i] == 0 && ideal_output[i] == 0){
      d=d+1; }
    
  }
  
  output = ((b+c)^2)/((a+b+c+d)^2);
  return(output)
}