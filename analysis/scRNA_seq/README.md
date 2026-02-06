Steps: 

Preprocessing 

1.	Create a data folder for the cellranger outputs for meninges, skull, and spine 
2.	run cellbender/soupX on all of them
3.	then demultiplex the ones that need demultiplexing (Soupercell) 
4.	then QC 
5.	Add harmonised metadata 

Integration 
6.	Integration scVI; covariates: organ, run_id (sequencing lane) or donor ID? 

Label Transfer: 

7.	Yolksac progenitors: https://developmental.cellatlas.io/yolk-sac
8.	Panfetal label transfer (Issac has the final annotations) 
Take the final yolksac object and transfer them to Ken’s data 


Notes: 

demultiplexing should be downstream of the cellbender/soupx
SOC needs to be run on cellbender barcodes (if using cellbender)
or cellranger-arc/cellranger barcodes (if using SOUPX)
