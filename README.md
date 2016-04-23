# io-touchpad

[![Build Status](https://travis-ci.org/0mp/io-touchpad.svg?branch=master)](https://travis-ci.org/0mp/io-touchpad) [![Code Health](https://landscape.io/github/0mp/io-touchpad/master/landscape.svg?style=plastic)](https://landscape.io/github/0mp/io-touchpad/master)

## Installation

    apt-get update
    apt-get install build-essential libatlas-dev libatlas3gf-base
    apt-get install python3-dev python3-setuptools python3-numpy python3-scipy python3-pip
    pip3 install scikit-learn
    make

## Usage

    cd app
    sudo ./app.py


## Description

 touchpadlib and app.py was implemented and linked.

 touchpadlib is library that contains functions responsible for taking the
   coordinates (and few other values of interest that will be used in later processing) of a touched point. It is
   modified evtest.

 app.py is basic application that uses touchpadlib - it is partitioned into two threads:
   one is continuously taking new event (info about touch) from touchpadlib and adds it to the queue (python structure
    used to synchronize threads)

   second is taking the oldest event from the queue and analyzes it. That means he either adds it to a collection of
    points representing symbol currently being drawn, or it determines that the symbol was finished. It does by checking
    if difference between signals is higher then 0.3s. (we're planning to add a possibility to recognize by a special
    "stop" signal (type of event)) It also cuts the stream of points if it is too long (in time) and sends what he has
    to symbol_interpreter. This is implemented so that if someone randomly draws shapes on touchpad without intention
    to send symbol, then they won't overflow memory. Currently symbol_interpreter can only count the number of events and
    write coordinates of first 10 received (it is used as a test for now).

### Usage

    make
    sudo python3 ./app.py

## Tests

### Installation

    apt-get update
    apt-get install python3-pytest

### Usage

    cd app/test
    py.test-3

## Tools

### matrixanalyser.py

#### Installation

    apt-get install python3-matplotlib
    cd app/tools
    make

#### Usage

    cd app
    sudo ./tools/matrixanalyser.py [--tolerance TOLERANCE]

All generated figures of the drawn symbols are stored inside
the `app/tools/data/matrixanalyser/figures` directory.
