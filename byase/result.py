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

from typing import Optional, NamedTuple, List
import io

import pandas as pd

from .db import DB, Field


class ResultRecord(NamedTuple):
    """Result record."""
    task_id: str
    success: bool
    error_msg: Optional[str]
    fragments_count: Optional[int]
    trace: Optional[pd.DataFrame]
    trace_stats: Optional[pd.DataFrame]


class ResultDB:
    """Inference result database."""

    _schema = [Field('id', 'S'), Field('success', 'I'), Field('error_msg', 'S'),
               Field('fragments_count', 'I'), Field('trace', 'B'), Field('trace_stats', 'B')]

    def __init__(self, db_path: str, initialize: bool = False, read_only: bool = True):
        self.db_path = db_path

        if initialize:
            assert not read_only

        self._db = DB(db_path, schema=ResultDB._schema if initialize else None, read_only=read_only)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._db.close()

    def get_all_task_ids(self):
        """Get all task IDs."""
        return self._db.get_all_item_ids()

    @staticmethod
    def _construct_record(item):
        """Construct record from database item."""
        task_id, success, error_msg, fragments_count, trace, stats = item
        if success == 0:
            success = False
            assert error_msg is not None
            assert fragments_count == -1
            assert trace is None
            assert stats is None
            fragments_count = None
        else:
            assert success == 1
            success = True
            assert error_msg is None

            trace_buf = io.BytesIO(trace)
            trace = pd.read_parquet(trace_buf)
            trace_buf.close()
            stats_buf = io.BytesIO(stats)
            stats = pd.read_parquet(stats_buf)
            stats_buf.close()
        return ResultRecord(task_id=task_id, success=success, error_msg=error_msg,
                            fragments_count=fragments_count, trace=trace, trace_stats=stats)

    def records_iterator(self):
        """Records iterator."""
        for item in self._db.item_iterator():
            yield self._construct_record(item)

    def records_db_fields_iterator(self):
        """Records database fields iterator."""
        for item in self._db.item_iterator():
            yield {self._schema[n].name: item[n] for n in range(len(self._schema))}

    def get_record(self, task_id):
        """Get record by task ID."""
        item = self._db.get_item(task_id)
        return self._construct_record(item)

    def store_record(self, record: ResultRecord):
        """Store result record."""
        if not record.success:
            assert record.error_msg is not None
            assert record.fragments_count is None
            assert record.trace is None
            assert record.trace_stats is None
            fragments_count, trace, stats = -1, None, None
        else:
            assert record.error_msg is None

            fragments_count = record.fragments_count

            trace_buf = io.BytesIO()
            record.trace.to_parquet(trace_buf, index=False)
            trace_buf.seek(0)
            trace = trace_buf.read()
            trace_buf.close()

            stats_buf = io.BytesIO()
            record.trace_stats.to_parquet(stats_buf, index=True)
            stats_buf.seek(0)
            stats = stats_buf.read()
            stats_buf.close()
        self._db.store_item((record.task_id, 1 if record.success else 0,
                             record.error_msg, fragments_count, trace, stats))

    def load_record(self, task_id: str):
        """Load result record."""
        item = self._db.get_item(task_id)
        return self._construct_record(item)

    def merge(self, paths: List[str]):
        """Merge other results."""
        for path in paths:
            with ResultDB(db_path=path) as other:
                for record in other.records_iterator():
                    self.store_record(record)
