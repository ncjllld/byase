# This file is part of BYASE.
#
# BYASE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BYASE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BYASE.  If not, see <https://www.gnu.org/licenses/>.
#
# Author: Lili Dong
#

import sys
from datetime import datetime


ERROR = 40
WARNING = 30
INFO = 20
DEBUG = 10


_LEVEL_NAME_MAP = {
    ERROR: 'ERROR',
    WARNING: 'WARNING',
    INFO: 'INFO',
    DEBUG: 'DEBUG'
}


class MessageCenter:
    """Message center.

    Attributes:
        _level: Logging level.
        _log_path: The path of log file.
    """

    def __init__(self, level=INFO, log_path=None):
        self._level = level
        self._log_path = log_path

    def _handle_log(self, log):
        """Log handler."""
        if self._log_path is None:
            print(log, file=sys.stderr)
        else:
            with open(self._log_path, 'a') as f:
                f.write(log + '\n')

    def _log(self, level, msg):
        """Log message based on given level."""
        if level < self._level:
            return None
        time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log = '{} | {} | {}'.format(time, _LEVEL_NAME_MAP[level], msg)
        self._handle_log(log)

    def log_error(self, msg):
        """Log error message."""
        self._log(ERROR, msg)

    def log_warning(self, msg):
        """Log warning message."""
        self._log(WARNING, msg)

    def log_info(self, msg):
        """Log info message."""
        self._log(INFO, msg)

    def log_debug(self, msg):
        """Log debug message."""
        self._log(DEBUG, msg)

    def handle_progress(self, msg, progress=None):
        """Handle progress"""
        if msg is not None:
            self.log_info(msg)
        else:
            assert progress is not None

    def handle_data(self, data):
        """Handle data."""
        pass
