% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_files.R
\name{path_append_and_check}
\alias{path_append_and_check}
\title{append a file to path. If longer than N characters (Windows filesystem max path length), try with md5 hashed filename}
\usage{
path_append_and_check(path, file)
}
\arguments{
\item{path}{a directory on this computer. eg; "C:/temp" on windows or "/home/user1" on unix}

\item{file}{filename that should be appended to the path}
}
\description{
note;
in contrast to using a timestamp or substring+index, this solution yields stable filenames. Should never fail if provided path/dir is reasonable
}
