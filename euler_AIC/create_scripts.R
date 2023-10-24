for(all_chi in 0:40){
  for (kk in 5:30){
      filecommand<-paste("\nk <- ",kk,"\nchi<-",all_chi,"\nsource(\"cluster_script.R\")\n",sep="")
      write(filecommand,file=paste("euler_scripts/cluster_euler_k", kk, "_chi",all_chi,".R",sep=""))
  }
}
