concatProfiles = function(p){
  if(is(p, "relSimDB")){
    result = as.vector(t(as.matrix(p$profiles)))
  }else if(is(p, "list") && is(p[[1]][[1]], "profile")){
    if(length(p) == 1){
      result = as.vector(t(do.call(rbind, p)))
    }else{
      result = as.vector(t(do.call(rbind, lapply(p, function(prof)do.call(rbind, prof)))))
    }
  }
  
  return(result)
}