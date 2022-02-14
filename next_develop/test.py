#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import os
import subprocess
import tempfile
import logging

import sys


VERSION = "1.0.1"


if sys.version_info[0] > 3:
    print("This script requires Python 2.7")
    exit(1)

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


    
import argparse
import json

subprocess.CalledProcessError(["which", "samtools"]).returncode

