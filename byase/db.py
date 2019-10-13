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

import os
from typing import List, NamedTuple
import struct


class DBError(Exception):
    """DB error."""
    def __init__(self, msg):
        super().__init__(msg)


class DBPathError(DBError):
    """DB path error."""
    def __init__(self, db_path, msg):
        super().__init__('[PATH: {}] {}'.format(db_path, msg))


class DBOperationError(DBError):
    """DB path error."""
    def __init__(self, db_path, msg):
        super().__init__('[DATABASE: {}] {}'.format(db_path, msg))


class Field(NamedTuple):
    name: str
    fmt: str


class ItemIdx(NamedTuple):
    id: str
    offset: int


_INT_FMT = '!q'
_INT_SIZE = struct.calcsize(_INT_FMT)


def _encode_int(value):
    """Encode int value."""
    return struct.pack(_INT_FMT, value)


def _read_int(file):
    """Read int from file."""
    buffer = file.read(_INT_SIZE)
    return struct.unpack(_INT_FMT, buffer)[0]


def _encode_field(field: Field, value):
    """Encode field value."""
    if field.fmt == 'I':
        assert isinstance(value, int)
        return _encode_int(value)

    if value is None:
        return _encode_int(0)

    if field.fmt == 'S':
        assert isinstance(value, str)
        value = value.encode()
    else:
        assert field.fmt == 'B'
        assert isinstance(value, bytes)
    return _encode_int(len(value)) + value


def _read_field(field: Field, file):
    """Read field from file."""
    if field.fmt == 'I':
        return _read_int(file)

    field_size = _read_int(file)
    if field_size == 0:
        return None

    buffer = file.read(field_size)
    if field.fmt == 'S':
        return buffer.decode()
    else:
        assert field.fmt == 'B'
        return buffer


def _encode_item(schema: List[Field], values):
    """Encode item."""
    assert len(schema) == len(values)
    code = b''
    for i, field in enumerate(schema):
        code += _encode_field(field, values[i])
    return _encode_int(len(code)) + code


def read_item(schema: List[Field], file):
    """Read item from file."""
    item_size = _read_int(file)
    start = file.tell()
    values = []
    for field in schema:
        value = _read_field(field, file)
        values.append(value)
    end = file.tell()
    assert item_size == end - start
    return values


class DB:
    """Database.

    Attributes:
        db_path: The path of database.
        read_only: If the database is read only.
        schema: The schema of database.
        idx_cache: The index cache.
    """

    _schema_info = [Field('field_name', 'S'), Field('field_fmt', 'S')]
    _idx_info = [Field('item_id', 'S'), Field('item_offset', 'I')]

    def __init__(self, db_path: str, schema: List[Field] = None, read_only: bool = True):
        self.db_path = db_path
        self.read_only = read_only

        if schema is not None:
            assert not self.read_only
            if os.path.exists(self.db_path):
                raise DBPathError(self.db_path, 'Database file already exists.')
            if os.path.exists(self.idx_path):
                raise DBPathError(self.idx_path, 'Database index file already exists.')
            self._store_schema(schema)

        self.schema = self._load_schema()

        if self.read_only:
            open_mode = 'rb'
        else:
            open_mode = 'ab+'

        self._db_file = open(self.db_path, open_mode)
        self._idx_file = open(self.idx_path, open_mode)

        if self.read_only:
            assert not self._db_file.writable()
            assert not self._idx_file.writable()
        else:
            assert self._db_file.writable()
            assert self._idx_file.writable()

        self.idx_cache = None

    def close(self):
        """Close."""
        self._db_file.close()
        self._idx_file.close()

    @property
    def idx_path(self):
        """The index path."""
        return self.db_path + '.idx'

    def _store_schema(self, schema: List[Field]):
        """Store field schema."""
        if self.read_only:
            raise DBOperationError(self.db_path, 'Database in read only mode, cannot store schema.')

        if 'id' not in [field.name for field in schema]:
            raise DBOperationError(self.db_path, 'Missing field "id" in schema.')
        for field in schema:
            if field.fmt not in ['I', 'S', 'B']:
                raise DBOperationError(self.db_path, 'Unknown fmt string "{}" in schema.'.format(field.fmt))

        schema_content = b''
        schema_size = 0
        for field in schema:
            content = _encode_item(DB._schema_info, (field.name, field.fmt))
            schema_content += content
            schema_size += len(content)
        with open(self.idx_path, 'wb') as file:
            file.write(_encode_int(schema_size))
            file.write(schema_content)

    def _load_schema(self):
        """Load field schema."""
        if not os.path.exists(self.idx_path):
            raise DBPathError(self.idx_path, 'Database index file not exists.')

        schema = []
        with open(self.idx_path, 'rb') as file:
            schema_size = _read_int(file)
            while file.tell() < _INT_SIZE + schema_size:
                name, fmt = read_item(DB._schema_info, file)
                schema.append(Field(name, fmt))
            assert file.tell() == _INT_SIZE + schema_size
        return schema

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _idx_iterator(self):
        """Index iterator."""
        f = self._idx_file
        f.seek(0, 2)
        file_size = f.tell()

        f.seek(0, 0)
        schema_size = _read_int(f)
        f.read(schema_size)
        while f.tell() < file_size:
            item = read_item(self._idx_info, f)
            item_id, item_offset = item
            yield ItemIdx(item_id, item_offset)

    def find_item_offset(self, item_id: str):
        """Find item offset."""
        if self.idx_cache is not None:
            return self.idx_cache[item_id] if item_id in self.idx_cache else None

        self.idx_cache = {}
        item_offset = None
        for idx in self._idx_iterator():
            self.idx_cache[idx.id] = idx.offset
            if item_id == idx.id:
                assert item_offset is None
                item_offset = idx.offset
        return item_offset

    def store_item(self, values):
        """Store item."""
        if self.read_only:
            raise DBOperationError(self.db_path, 'Database in read only mode, cannot store item.')

        item_id_idx = [field.name for field in self.schema].index('id')
        item_id = values[item_id_idx]

        if self.find_item_offset(item_id) is not None:
            raise DBOperationError(self.db_path, 'Item with ID "{}" already exists in database.'.format(item_id))

        self._db_file.seek(0, 2)
        item_offset = self._db_file.tell()

        # Store item.
        content = _encode_item(self.schema, values)
        self._db_file.write(content)
        # Store item index.
        content = _encode_item(self._idx_info, (item_id, item_offset))
        self._idx_file.write(content)

        # Update index cache.
        if self.idx_cache is not None:
            self.idx_cache[item_id] = item_offset

    def get_item_by_offset(self, item_offset: int):
        """Get item by offset."""
        self._db_file.seek(item_offset)
        values = read_item(self.schema, self._db_file)
        return values

    def get_item(self, item_id: str):
        """Get item by ID."""
        item_offset = self.find_item_offset(item_id)
        if item_offset is None:
            raise DBOperationError(self.db_path, 'Item with ID "{}" not exists in database.'.format(item_id))
        return self.get_item_by_offset(item_offset)

    def item_iterator(self):
        """Item iterator."""
        for idx in self._idx_iterator():
            yield self.get_item_by_offset(idx.offset)

    def get_all_item_ids(self):
        """Get all item IDs."""
        return [idx.id for idx in self._idx_iterator()]

    def preview(self, n: int):
        """Preview items."""
        print([field.name for field in self.schema])
        m = -1
        for k, item in enumerate(self.item_iterator()):
            m += 1
            if m >= n:
                break
            for i, v in enumerate(item):
                if isinstance(v, bytes):
                    item[i] = v[:min(len(v), 15)]
            print(k, item)
