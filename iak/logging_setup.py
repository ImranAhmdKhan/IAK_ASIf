from __future__ import annotations

import logging
import queue
import sys
from typing import Optional


class _QueueHandler(logging.Handler):
    def __init__(self, q: queue.Queue):
        super().__init__()
        self._q = q

    def emit(self, record: logging.LogRecord):
        self._q.put(self.format(record))


class SafeStreamHandler(logging.StreamHandler):
    def emit(self, record):
        try:
            msg = self.format(record)
            self.stream.write(msg + self.terminator)
            self.flush()
        except UnicodeEncodeError:
            msg = self.format(record)
            safe_msg = msg.encode(sys.stdout.encoding or "ascii", "replace").decode(sys.stdout.encoding or "ascii")
            try:
                self.stream.write(safe_msg + self.terminator)
                self.flush()
            except Exception:
                pass
        except Exception:
            self.handleError(record)


def setup_logging(verbose: bool = False, gui_queue: Optional[queue.Queue] = None):
    level = logging.DEBUG if verbose else logging.INFO
    logger = logging.getLogger("IAK")
    logger.setLevel(level)
    logger.handlers.clear()
    ch = SafeStreamHandler(sys.stdout)
    formatter = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    if gui_queue is not None:
        qh = _QueueHandler(gui_queue)
        qh.setFormatter(formatter)
        logger.addHandler(qh)
    return logger
