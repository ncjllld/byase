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

from .annotation import AnnotationDB
from .inference import InferenceTool
from .result import ResultDB
from .task.result import TaskResult
from .task.plot import TaskPlot, TaskPlotType, TaskPlotError


_PLOT_DIR_NAME = 'plots'


class PlotError(Exception):
    """Plot error."""
    def __init__(self, msg):
        super().__init__(msg)


class PlotPathError(PlotError):
    """Plot path error."""
    def __init__(self, path: str, msg):
        super().__init__('[PLOT ERROR] [PATH: {}] {}'.format(path, msg))


def plot_task(args):
    """Plot task."""
    result_dir = args['result_dir']
    task_id = args['task_id']
    mc = args['mc']

    inference_tool = InferenceTool(result_dir, n_process=1, param=None, mc=mc)

    with AnnotationDB(inference_tool.param.anno_path) as anno_db:
        task = anno_db.get_task(task_id)

    with ResultDB(inference_tool.result_path) as result_db:
        record = result_db.load_record(task_id)

    if not record.success:
        raise TaskPlotError(task_id, 'The inference for the task wss failed.')

    task_result = TaskResult(task, record.trace)

    plot_dir = os.path.join(result_dir, _PLOT_DIR_NAME)
    task_plot_dir = os.path.join(plot_dir, task_id)
    for dir_path in [plot_dir, task_plot_dir]:
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)

    plot_type = TaskPlotType.ASE

    task_plot = TaskPlot(task, inference_tool.bam_param,
                         task_result.trace, record.trace_stats,
                         task_plot_dir, plot_type, mc)
    html_path = task_plot.plot()
    return html_path
