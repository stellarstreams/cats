#!/bin/zsh
rsync -av -e ssh --exclude='poly*' --exclude='mast-*' rusty:"~/projects/cats/data/" ~/projects/cats/data/
