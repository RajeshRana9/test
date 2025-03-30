#!/bin/bash
set -e  # Exit on error
sudo apt-get update
sudo apt-get install -y \
    libxml2-dev \
    libxslt1-dev \
    zlib1g-dev \
    libgl1-mesa-glx
