#!/bin/bash

framerate=30
plotfile=plot.plt

while [ "$1" != "" ]; do
    case $1 in
        -f ) shift
            framerate=$1
            ;;
        -p ) shift
            plotfile=$1
            ;;
    esac
    shift
done

gnuplot "$plotfile" || exit 1
(pushd frames; ffmpeg -y -framerate "$framerate" -i frame_%05d.png video.mp4)
#mpv frames/video.mp4
