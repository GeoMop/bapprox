B-Spline Approximator
=====================

[![Build Status](https://travis-ci.org/GeoMop/bapprox.svg?branch=master)](https://travis-ci.org/GeoMop/bapprox)
[![Code Health](https://landscape.io/github/GeoMop/bapprox/master/landscape.svg?style=flat)](https://landscape.io/github/GeoMop/bapprox/master)
[![Code Climate](https://codeclimate.com/github/GeoMop/bapprox/badges/gpa.svg)](https://codeclimate.com/github/GeoMop/bapprox)
[![Test Coverage](https://codeclimate.com/github/GeoMop/bapprox/badges/coverage.svg)](https://codeclimate.com/github/GeoMop/bapprox/coverage)

Python tool for approximate 3d points using Bspline surface.

![Bapprox screenshot](/data/screenshot.png "Bapprox screenshot")

Requirements
------------

* Python 2.7
* Numpy 1.11
* SciPy 0.18
* GnuPlot 5.0
* Python OCC 0.17

How to Run
----------

    python ./src/main.py -f ./conf.yaml

How to Test
-----------

    PYTHONPATH=$PWD/src py.test tests

