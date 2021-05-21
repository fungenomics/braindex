# Script to populate the data directory for this app.
#   - Reads data.json (metadata on each required file)
#   - Creates necessary subdirectories in destination
#   - Identifies which files need to be copied over from hydra
#   - Prints a command to copy the files (using rsync)
#
# Usage:
# 1. Run the script from within the destination directory with:
#    $ Rscript populate_data_small_vers.R <first.last>
#    where <first.last> is your hydra username (the section 
#    before @hydra.ladydavis.ca)
#
# 2. Copy and execute the rsync commands produced as output

library(dplyr)
library(tidyr)
library(magrittr)
library(rjson)
library(purrr)
library(glue)

# Store hydra username - command line argument
args = commandArgs()

# Check if username was provided
if (length(args) == 0) {
  stop("Please supply your username")
} 

path_to_projects <- "/home/kleinman/"

# Read in data.json, returning a list with one element per folder
app_data <- rjson::fromJSON(file = "data.json")

# Get the names of the directories
directories <- names(app_data)

# Create directories
walk(directories, dir.create, showWarnings = FALSE)

# For each list,
#     look for the files which are not produced by a script, and copy those over
process_file <- function(file, dir_name) {
  
  # Copy file locally only if it is not produced by a script
  if (is.null(file$script)) {
    
    # If a small version of the file does not exist, use the path to full file
    if (is.null(file$small_ver_path)) {
      
      src  <- file.path(path_to_projects, file$path)
      
    } else {
      
      # Use the path to the small version for local use
      src <- file.path(path_to_projects, file$small_ver_path)
      
    }
    
    # Produce destination path based on current directory
    dest <- file.path(dir_name, file$file)
    
    # Produce and echo rsync command including hydra username
    cmd <- glue("rsync {args[1]}@hydra.ladydavis.ca:{src} {dest}")
    message(cmd)
    
  } else {
    
    NULL
    
  }
  
}

# Loop over the directories
iwalk(app_data, function(dir_data, dir_name) {
  
  # Loop over files, and get the copy commands, if it doesn't have a script property
  walk(dir_data, ~ process_file(.x, dir_name = dir_name))
  
})