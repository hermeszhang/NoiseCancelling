#!/bin/sh

export LD_LIBRARY_PATH=/opt/altera/11.1/quartus/linux

if [ -x /opt/grmon-eval/linux/grmon-eval ]; then
    /opt/grmon-eval/linux/grmon-eval -altjtag -u
fi

unset LD_LIBRARY_PATH
