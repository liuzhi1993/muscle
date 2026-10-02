# Install and switch between the local MUSCLE package versions.

MUSCLE_TARBALLS <- c(
  `0.1.1` = "/Users/liuzhi1993/Desktop/MUSCLE_Paper_Correction/muscle_code/muscle/muscle_0.1.1.tar.gz",
  `0.1.2` = "/Users/liuzhi1993/Desktop/MUSCLE_Paper_Correction/muscle_0.1.2.tar.gz",
  `0.1.3` = "/Users/liuzhi1993/Desktop/MUSCLE_Paper_Correction/muscle-github/muscle_0.1.3.tar.gz"
)

MUSCLE_VERSION_LIB <- file.path(
  "/Users/liuzhi1993/Desktop/MUSCLE_Paper_Correction",
  ".muscle_versions"
)

available_muscle_versions <- function() {
  names(MUSCLE_TARBALLS)
}

install_muscle_version <- function(version, force = FALSE) {
  version <- as.character(version)
  if (!version %in% names(MUSCLE_TARBALLS)) {
    stop("version must be one of: ", paste(names(MUSCLE_TARBALLS), collapse = ", "))
  }

  tarball <- MUSCLE_TARBALLS[[version]]
  if (!file.exists(tarball)) {
    stop("Package tarball does not exist: ", tarball)
  }

  lib <- file.path(MUSCLE_VERSION_LIB, version)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  installed_version <- tryCatch(
    as.character(utils::packageVersion("muscle", lib.loc = lib)),
    error = function(e) NA_character_
  )

  if (isTRUE(!force) && identical(installed_version, version)) {
    return(invisible(lib))
  }

  utils::install.packages(
    tarball,
    repos = NULL,
    type = "source",
    lib = lib
  )

  installed_version <- as.character(utils::packageVersion("muscle", lib.loc = lib))
  if (!identical(installed_version, version)) {
    stop("Installed version ", installed_version, " instead of ", version)
  }
  invisible(lib)
}

switch_muscle_version <- function(version, force_install = FALSE) {
  version <- as.character(version)
  lib <- install_muscle_version(version, force = force_install)

  if ("package:muscle" %in% search()) {
    detach("package:muscle", unload = TRUE, character.only = TRUE)
  }
  if ("muscle" %in% loadedNamespaces()) {
    try(unloadNamespace("muscle"), silent = TRUE)
  }

  .libPaths(unique(c(lib, .libPaths())))
  suppressPackageStartupMessages(
    library("muscle", character.only = TRUE, lib.loc = lib)
  )

  active_version <- as.character(utils::packageVersion("muscle"))
  if (!identical(active_version, version)) {
    stop("Active version is ", active_version, " instead of ", version)
  }
  message("Active MUSCLE version: ", active_version)
  invisible(active_version)
}

