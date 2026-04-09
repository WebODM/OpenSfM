## OpenSfM

> Meta has not been actively developing [OpenSfM](https://github.com/mapillary/opensfm). We're hoping the fantastic work of OpenSfM can continue to live on here at [WebODM](https://github.com/WebODM). We're ready to support development efforts.

OpenSfM ![Docker workflow](https://github.com/WebODM/opensfm/workflows/Docker%20CI/badge.svg)


=======

## Overview
OpenSfM is a Structure from Motion library written in Python. The library serves as a processing pipeline for reconstructing camera poses and 3D scenes from multiple images. It consists of basic modules for Structure from Motion (feature detection/matching, minimal solvers) with a focus on building a robust and scalable reconstruction pipeline. It also integrates external sensor (e.g. GPS, accelerometer) measurements for geographical alignment and robustness. A JavaScript viewer is provided to preview the models and debug the pipeline.

<p align="center">
  <img src="https://opensfm.org/docs/_images/berlin_viewer.jpg" />
</p>

Checkout this [blog post with more demos](http://blog.mapillary.com/update/2014/12/15/sfm-preview.html)


## Getting Started

* [Building the library][]
* [Running a reconstruction][]
* [Documentation][]


[Building the library]: https://opensfm.org/docs/building.html (OpenSfM building instructions)
[Running a reconstruction]: https://opensfm.org/docs/using.html (OpenSfM usage)
[Documentation]: https://opensfm.org/docs/ (OpenSfM documentation)

## License

The original OpenSfM code is BSD-style licensed.

The current code in this repository is AGPLv3 licensed, as found in the LICENSE file.
