saveModelObj <- function(model_obj, save_name = "masc.modelobj.rds", save_dir = NULL) {
  # If directory is unspecified, use current working directory
  if(is.null(save_dir)) {
    message("save_model_dir is unspecified, saving model objects to current directory")
    save_dir <- getwd()
  }
  # Try to save file unless it already exists
  if(file.exists(save_name) == FALSE) {
    saveRDS(model_obj, file = save_name)
    message(paste("Models saved to", file.path(save_dir, save_name)))
  } else {
    warning(paste(save_name, "already exists in directory, did not overwrite"))
  }
}
