## Test environments
* local OS X install, R 3.6.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was no NOTE on local OS X.

There was 1 NOTE on win-builder:

* Possibly mis-spelled words in DESCRIPTION:
  Yu (7:76)
  al (7:82, 7:332)
  et (7:79, 7:329)
  untruncated (7:282)

  "Yu et al" and "Lin et al" refer to the papers this implementation is based on. The word "untruncated" is the antonym of "truncated". These are therefore not mis-spelled words.
