#!/usr/bin/env python
# Copyright (C) 2010 Vinay Sajip. All Rights Reserved.
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted,
# provided that the above copyright notice appear in all copies and that
# both that copyright notice and this permission notice appear in
# supporting documentation, and that the name of Vinay Sajip
# not be used in advertising or publicity pertaining to distribution
# of the software without specific, written prior permission.
# VINAY SAJIP DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
# ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL
# VINAY SAJIP BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
# ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
# IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
# OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
"""
An example script showing how to use logging with multiprocessing.

The basic strategy is to set up a listener process which can have any logging
configuration you want - in this example, writing to rotated log files. Because
only the listener process writes to the log files, you don't have file
corruption caused by multiple processes trying to write to the file.

The listener process is initialised with a queue, and waits for logging events
(LogRecords) to appear in the queue. When they do, they are processed according
to whatever logging configuration is in effect for the listener process.

Other processes can delegate all logging to the listener process. They can have
a much simpler logging configuration: just one handler, a QueueHandler, needs
to be added to the root logger. Other loggers in the configuration can be set
up with levels and filters to achieve the logging verbosity you need.

A QueueHandler processes events by sending them to the multiprocessing queue
that it's initialised with.

In this demo, there are some worker processes which just log some test messages
and then exit.

This script was tested on Ubuntu Jaunty and Windows 7.

Copyright (C) 2010 Vinay Sajip. All Rights Reserved.
"""
# You'll need these imports in your own code
import logging
import logging.handlers
from logging.config import fileConfig
from multiprocessing import Queue


class QueueHandler(logging.Handler):
    """
    This is a logging handler which sends events to a multiprocessing queue.

    The plan is to add it to Python 3.2, but this can be copy pasted into
    user code for use with earlier Python versions.
    """

    def __init__(self, queue: Queue):
        """
        Initialise an instance, using the passed queue.
        """
        logging.Handler.__init__(self)
        self.queue = queue

    def emit(self, record):
        """
        Emit a record.

        Writes the LogRecord to the queue.
        """
        # noinspection PyBroadException
        try:
            ei = record.exc_info
            if ei:
                _ = self.format(record)  # just to get traceback text into record.exc_text
                record.exc_info = None  # not needed any more
            self.queue.put_nowait(record)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


def log_listener(queue: Queue, config_file: str):
    # logging config
    fileConfig(config_file)

    while True:
        # noinspection PyBroadException
        try:
            record = queue.get()
            if record is None:  # We send this as a sentinel to tell the listener to quit.
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)  # No level or filter logic applied - just do it!
        except Exception:
            import sys
            import traceback
            print('Whoops! Problem:', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)


def log_listener_no_config(queue: Queue):
    while True:
        # noinspection PyBroadException
        try:
            record = queue.get()
            if record is None:  # We send this as a sentinel to tell the listener to quit.
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)  # No level or filter logic applied - just do it!
        except Exception:
            import sys
            import traceback
            print('Whoops! Problem:', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
