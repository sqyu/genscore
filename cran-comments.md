## Resubmission
Fixed memory issues.

## Test environments
* local OS X install, R 3.6.3
* win-builder (devel and release)
* R-hub builder:
  - Debian Linux, R-devel, clang, ISO-8859-15 locale
  - Debian Linux, R-devel, GCC
  - Debian Linux, R-devel, GCC ASAN/UBSAN
  - Debian Linux, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
  - Fedora Linux, R-devel, GCC
  - Ubuntu Linux 16.04 LTS, R-devel, GCC
  - Ubuntu Linux 16.04 LTS, R-release, GCC


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE in some environments, and none in the others:

* Possibly mis-spelled words in DESCRIPTION:
  al (7:82, 7:373)
  et (7:79, 7:370)
  untruncated (7:323)
  Yu (7:76)

