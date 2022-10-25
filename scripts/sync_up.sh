#!/bin/zsh
rsync -av -e ssh --exclude='poly*' --exclude='mast-*' ~/projects/cats/data/ rusty:"~/projects/cats/data/"
