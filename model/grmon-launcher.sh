#!/bin/sh

export LD_LIBRARY_PATH=/opt/altera/11.1/quartus/linux

if [ -x /opt/grmon-eval-2.0.36/linux/bin/grmon ]; then
    /opt/grmon-eval-2.0.36/linux/bin/grmon -altjtag -u
fi

unset LD_LIBRARY_PATH
