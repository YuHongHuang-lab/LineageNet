execute_py = function(python_path=NULL,
                      script_path=NULL,
                      ...){
  args = c(script_path,
           unname(unlist(list(...))))
  system2(python_path, args = args, wait = T)
}
