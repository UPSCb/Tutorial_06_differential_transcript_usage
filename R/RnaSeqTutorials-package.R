# RnaSeqTutorial06 package
#
# * To find the imported packages, in the terminal
#
# ---
# cd inst
# grep "library(" */*/*.Rmd | sed -e 's:.*library::g' | tr -d '()' | sort | uniq
# ---
#
# * To build the DESCRIPTION Imports string
#
# ---
# library(here)
# pkgs <- c("DRIMSeq","EnsDb.Hsapiens.v86","ggvenn","here","learnr","org.Hs.eg.db","stageR","tidyverse","tximport")
# write(paste0("    ",pkgs," (>= ",unlist(installed.packages()[pkgs,"Version"],use.names=FALSE),"),"),
#       file="Imports.tmp")
# ---
#
#' @title RnaSeqTutorials package
#' @section Tutorials:
#' This is the sixth in a series of tutorials
#' \itemize{
#' \item\code{06_differential_transcript_usage} a tutorial on differential transcript usage
#' }
#'
#' @name RnaSeqTutorial06 package
#' @rdname RnaSeqTutorial06-package
#' @author Nicolas Delhomme [aut,cre]
#' @keywords package
#' @description A simple description of the RnaSeqTutorial06 package
#' @seealso The vignette
#' @examples
#' 	\dontrun{
#' 	learnr::run_tutorial("06_differential_transcript_usage", package = "RnaSeqTutorial06")
#' 	}
#' @keywords internal
"_PACKAGE"
#'
NULL
