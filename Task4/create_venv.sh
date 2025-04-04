#!/bin/bash
/usr/bin/python3 -m venv .venv &&
        source .venv/bin/activate &&
        pip install --upgrade pip setuptools wheel &&
        pip3 install pandas numpy matplotlib
