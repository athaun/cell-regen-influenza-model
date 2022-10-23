#!/usr/bin/env bash
clear
echo "Running test"
screen -S viralTransGPU bash -c "tmux new-session \; \
                          send-keys './multirun.sh' C-m \; \
                          split-window -h \; \
                          send-keys 'htop' C-m \; \
                          split-window -v \;"


                          
