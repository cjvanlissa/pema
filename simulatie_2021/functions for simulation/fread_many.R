fread_many = function(files,header=F,...){
  if(length(files)==0) return()
  if(typeof(files)!='character') return()
  files = files[file.exists(files)]
  if(length(files)==0) return()
  tmp = tempfile(fileext = ".csv")
  if(header==T){
    system(paste0('head -n1 ',files[1],' > ',tmp))
    system(paste0("xargs awk 'FNR>1' >> ",tmp),input=files)
  } else {
    system(paste0("xargs awk '1' > ",tmp),input=files)
  }
  DT = fread(file=tmp,header=header,...)
  file.remove(tmp)
  DT
}
