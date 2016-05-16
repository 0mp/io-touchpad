# io-touchpad

[![Build Status](https://travis-ci.org/0mp/io-touchpad.svg?branch=master)](https://travis-ci.org/0mp/io-touchpad)
[![Code Health](https://landscape.io/github/0mp/io-touchpad/master/landscape.svg?style=flat)](https://landscape.io/github/0mp/io-touchpad/master)
[![Code Issues](https://www.quantifiedcode.com/api/v1/project/538d6bb306774bd7ae52b8c4dbdd0854/badge.svg)](https://www.quantifiedcode.com/app/project/538d6bb306774bd7ae52b8c4dbdd0854)
[![codecov](https://codecov.io/gh/0mp/io-touchpad/branch/master/graph/badge.svg)](https://codecov.io/gh/0mp/io-touchpad)


## Application

### Installation

    apt-get update
    apt-get install build-essential libatlas-dev libatlas3gf-base
    apt-get install python3-dev python3-setuptools python3-numpy python3-scipy python3-pip
    pip3 install scikit-learn
    make

### Usage

1. Go to the `app` directory.

        cd app

2. Undertake a learning session with a trainig size of SIZE for a symbol named SYMBOL:

        sudo ./app.py learn SIZE SYMBOL

3. Use the app. (Hopefully, it will recognise your symbol.)

        sudo ./app run user

If you want to use predefined, hardcoded symbols run `./app --help` for more information.



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
    sudo ./tools/matrixanalyser.py [--help] [--tolerance TOLERANCE] [--show]

More information on the usege is avaialabe if you run `./tools/matrixanalyser.py --help`.


All generated figures of the drawn symbols are stored inside
the `app/tools/data/matrixanalyser/figures` directory.
