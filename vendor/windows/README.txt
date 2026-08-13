xfun 0.55 Windows compatibility binary

Purpose:
The current CRAN R 4.4 Windows binary repository can install knitr 1.49 with
xfun 0.57. knitr 1.49 imports xfun::attr, which xfun 0.57 no longer exports.
This checked xfun 0.55 binary restores the compatible API for R 4.4 users.

Upstream source:
https://cran.r-project.org/src/contrib/Archive/xfun/xfun_0.55.tar.gz

Source SHA-256:
398fc5136d3b8ca8d09bd5987e8f10421dec77f0e1175704a9f5f2d1ceb5d36e

Binary SHA-256:
6c284c66db4d88cc4306f68a6ccce0a8b07e9d0ecc5bd1c6d1093b13f81c9ad4

Build environment:
R 4.4.2 for Windows with Rtools44, using R CMD INSTALL --build.

License:
MIT + file LICENSE. The installed binary retains the package DESCRIPTION,
LICENSE, help, and metadata supplied by the upstream source.
