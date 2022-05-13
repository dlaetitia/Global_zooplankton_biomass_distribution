#
# Common settings and functions
#
# (c) 2022 L Drago, T Panaïotis, J-O Irisson, GNU General Public License v3

suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
library("ecotaxar")
#remotes::install_github("jiho/chroma", force = TRUE)
library("chroma")
library("ggplot2")


# source for large data that does not fit in the project
data_dir <- config$data_dir

# number of cores for parallel computation
n_cores <- unique(config$n_cores)

# create two functions to save figures on their own through savefig() or in a pdf document through openfig()
savefig <- function(name, width=18, height=20, units="cm", ...) {
  ggsave(filename=paste0(figure_dir, name), width=width, height=height, units=units, ...)
}
openfig <- function(name, width=18, height=20, ...) {
  pdf(file=paste0(figure_dir, name), width=width/2.54, height=height/2.54, ...)
}

# load coast data
coast <- read_csv(paste0("Data/", config$element, "/coast_world.csv", col_types=cols()))