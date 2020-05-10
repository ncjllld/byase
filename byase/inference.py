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
from multiprocessing import Process, Queue
from typing import NamedTuple, List, Optional
import shutil

from .message import MessageCenter
from .task import BAMParam
from .annotation import AnnotationDB
from .result import ResultRecord, ResultDB


_PARAM_FILENAME = 'param.tab'
_RESULT_FILENAME = 'result.db'
_TMP_DIR_NAME = 'tmp'
_SIGNAL_END = -1


class InferenceError(Exception):
    """Inference error."""
    def __init__(self, msg):
        super().__init__(msg)


class InferencePathError(InferenceError):
    """"Path error."""
    def __init__(self, path, msg):
        super().__init__('[PATH: {}] {}'.format(path, msg))


class InferenceParseError(InferenceError):
    """"Parsing error."""
    def __init__(self, path, msg):
        super().__init__('[PARSE: {}] {}'.format(path, msg))


def _inference_task(task_id: str, anno_path: str, bam_param: BAMParam, mcmc_samples: int, tune_samples: int, q: Queue):
    """Inference for task."""
    os.environ['OMP_NUM_THREADS'] = '1'

    from .task.inference import TaskInference, TaskInferenceError
    from .task.result import TaskResult

    mc = MessageCenter()

    with AnnotationDB(anno_path) as anno_db:
        task = anno_db.get_task(task_id)

    error = None
    try:
        task_inference = TaskInference(task, bam_param, mcmc_samples, tune_samples, mc)
        fragments_count, trace = task_inference.run()
        task_result = TaskResult(task, trace)
        stats = task_result.stats()
    except Exception as e:
        fragments_count = None
        trace = None
        stats = None
        error = '{}'.format(e) if isinstance(e, TaskInferenceError) else '[UNKNOWN ERROR] {}'.format(e)
        mc.log_error('[TASK: {}] {}'.format(task_id, error))
    success = error is None

    record = ResultRecord(task_id=task_id, success=success, error_msg=error,
                          fragments_count=fragments_count, trace=trace, trace_stats=stats)

    q.put(record)


def _handle_tasks(num: int, anno_path: str, bam_param: BAMParam, mcmc_samples: int, tune_samples: int,
                  theano_cache_dir: str, in_queue: Queue, out_queue: Queue, mc: MessageCenter):
    """Handle tasks."""
    theano_cache_path = os.path.join(theano_cache_dir, str(num))
    os.environ['THEANO_FLAGS'] = 'base_compiledir={}'.format(theano_cache_path)

    q = Queue()

    n = 0
    while True:
        item = in_queue.get()
        if isinstance(item, int) and item == _SIGNAL_END:
            break
        assert isinstance(item, str)

        task_id = item
        mc.handle_progress('[PROCESS {}] Inference for task: [{}]...'.format(num, task_id))
        mc.handle_data(('process_start', task_id))

        proc = Process(target=_inference_task, args=(task_id, anno_path, bam_param, mcmc_samples, tune_samples, q))
        proc.start()
        out_queue.put(q.get())
        proc.join()

        n += 1
        if n % 100 == 0:
            if os.path.exists(theano_cache_path):
                mc.log_debug('[PROCESS {}] Clear theano compiling cache...'.format(num))
                shutil.rmtree(theano_cache_path)
    out_queue.put(_SIGNAL_END)


def _dispatch_tasks(task_ids: List[str], handlers_count: int, in_queue: Queue):
    """Dispatch tasks."""
    for task_id in task_ids:
        in_queue.put(task_id)
    for _ in range(handlers_count):
        in_queue.put(_SIGNAL_END)


def _store_result_record(result_path: str, record: ResultRecord):
    """Store result record."""
    with ResultDB(result_path, read_only=False) as result_db:
        result_db.store_record(record)


class InferenceParam(NamedTuple):
    anno_path: str
    bam_param: BAMParam
    mcmc_samples_count: int
    tune_samples_count: int


class InferenceTool:
    """Inference tool.

    Attributes:
        n_process: The count of processes.
        param: The inference param.
        out_dir: The path of output directory.
        mc: Message center.
    """

    def __init__(self, out_dir: str, n_process: int, param: Optional[InferenceParam], mc: MessageCenter):
        self.out_dir = out_dir
        self.n_process = n_process
        self.mc = mc

        if not os.path.isdir(self.out_dir):
            raise InferencePathError(out_dir, 'Invalid directory.')

        fresh_start = param is not None

        if fresh_start:
            self.param = param
            if os.path.exists(self.result_path):
                raise InferencePathError(self.result_path, 'The result file already exists.')
        else:
            self.param = self.load_param()

        # Check paths in params.
        if not os.path.isdir(self.anno_path):
            raise InferencePathError(self.anno_path, 'The annotation directory is invalid.')
        for path in self.bam_param.paths:
            if not os.path.exists(path):
                raise InferencePathError(path, 'The BAM file not exists.')

        self.mc.log_debug('n_process: {}'.format(self.n_process))
        self.mc.log_debug('mcmc_samples: {}'.format(self.param.mcmc_samples_count))
        self.mc.log_debug('tune_samples: {}'.format(self.param.tune_samples_count))
        self.mc.log_debug('out_dir: {}'.format(self.out_dir))
        self.mc.log_debug('anno_path: {}'.format(self.anno_path))
        self.mc.log_debug('bam_paths: {}'.format(self.bam_param.paths))
        self.mc.log_debug('read_len: {}'.format(self.bam_param.read_len))
        self.mc.log_debug('paired_end: {}'.format(self.bam_param.paired_end))
        self.mc.log_debug('insert_size_mean: {}'.format(self.bam_param.insert_size_mean))
        self.mc.log_debug('insert_size_std: {}'.format(self.bam_param.insert_size_std))

        if fresh_start:
            # Store params.
            self._store_param()

        if not os.path.exists(self.result_path):
            with ResultDB(self.result_path, initialize=True, read_only=False):
                pass

    @property
    def anno_path(self):
        """The path of annotation."""
        return self.param.anno_path

    @property
    def bam_param(self):
        """BAM params."""
        return self.param.bam_param

    @property
    def param_path(self):
        """The path of params."""
        return os.path.join(self.out_dir, _PARAM_FILENAME)

    @property
    def result_path(self):
        """The path of results."""
        return os.path.join(self.out_dir, _RESULT_FILENAME)

    @property
    def theano_cache_dir(self):
        """The path of Theano cache directory."""
        return os.path.join(os.path.abspath(self.out_dir), 'theano_cache')

    @property
    def tmp_dir(self):
        """The path of tmp directory."""
        return os.path.join(self.out_dir, _TMP_DIR_NAME)

    @property
    def tmp_result_path(self):
        """The path of tmp result path."""
        return os.path.join(self.tmp_dir, _RESULT_FILENAME)

    def _store_param(self):
        """Store params."""
        if os.path.exists(self.param_path):
            raise InferencePathError(self.param_path, 'The inference parameters file already exists.')

        bp = self.bam_param
        data = [('Annotation', os.path.relpath(self.anno_path, self.out_dir)),
                ('BAM', ';'.join([os.path.relpath(path, self.out_dir) for path in bp.paths])),
                ('MCMC-samples', self.param.mcmc_samples_count),
                ('Tune-samples', self.param.tune_samples_count),
                ('Read-len', bp.read_len)]
        if bp.paired_end:
            data.append(('Paired-end', '{};{}'.format(bp.insert_size_mean, bp.insert_size_std)))
        with open(self.param_path, 'w') as f:
            for k, v in data:
                f.write('{}\t{}\n'.format(k, v))

    def load_param(self):
        """Load params."""
        if not os.path.exists(self.param_path):
            raise InferencePathError(self.param_path, 'The inference parameters file not exists.')

        anno_path, bam_paths, read_len = None, None, None
        mcmc_samples_count, tune_samples_count = None, None
        paired_end, insert_size_mean, insert_size_std = False, None, None

        with open(self.param_path) as f:
            for line in f:
                line = line.strip('\n')
                k, v = line.split('\t')
                if k == 'Annotation':
                    anno_path = os.path.join(self.out_dir, v)
                elif k == 'BAM':
                    bam_paths = [os.path.join(self.out_dir, path) for path in v.split(';')]
                elif k == 'MCMC-samples':
                    mcmc_samples_count = int(v)
                elif k == 'Tune-samples':
                    tune_samples_count = int(v)
                elif k == 'Read-len':
                    read_len = int(v)
                elif k == 'Paired-end':
                    paired_end = True
                    m, n = v.split(';')
                    insert_size_mean, insert_size_std = float(m), float(n)

        if (anno_path is None) or (bam_paths is None) or (read_len is None):
            raise InferenceParseError(self.param_path, 'The file contains incomplete inference parameters.')

        bam_param = BAMParam(bam_paths=bam_paths, read_len=read_len, paired_end=paired_end,
                             insert_size_mean=insert_size_mean, insert_size_std=insert_size_std)

        return InferenceParam(anno_path=anno_path, bam_param=bam_param,
                              mcmc_samples_count=mcmc_samples_count, tune_samples_count=tune_samples_count)

    def _store_tmp_result_records(self, records_queue: Queue):
        """Store tmp result records."""
        n = 0
        finished_count = 0
        while True:
            item = records_queue.get()
            if isinstance(item, int) and item == _SIGNAL_END:
                finished_count += 1
                if finished_count == self.n_process:
                    break
                continue
            assert isinstance(item, ResultRecord)

            proc = Process(target=_store_result_record, args=(self.tmp_result_path, item))
            proc.start()
            proc.join()

            n += 1
            self.mc.handle_data(('process_end', item.task_id))
            self.mc.handle_progress('Task: [{}] finished. ({} tasks processed)'.format(item.task_id, n), progress=n)

    def _merge_results(self):
        """Merge results."""
        self.mc.handle_progress('Merging results...')

        with ResultDB(self.result_path, read_only=False) as result_db:
            result_db.merge([self.tmp_result_path])

    def _extract_task_ids(self, target_count: Optional[int]):
        """Extract task IDs to be processed."""
        with AnnotationDB(self.anno_path) as anno_db:
            all_task_ids = anno_db.get_all_task_ids()
        with ResultDB(self.result_path) as result_db:
            finished_task_ids = set(result_db.get_all_task_ids())

        n = 0
        task_ids = []
        for task_id in all_task_ids:
            if task_id not in finished_task_ids:
                task_ids.append(task_id)
                n += 1
                if (target_count is not None) and (n == target_count):
                    break
        return task_ids

    def run(self, target_count: Optional[int], target_task_ids: Optional[List[str]] = None):
        """Run inference."""
        if os.path.exists(self.tmp_dir):
            self.mc.log_warning('The last run may not be completed, the tmp directory is not removed.')
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)
        with ResultDB(self.tmp_result_path, initialize=True, read_only=False):
            pass

        if target_task_ids is None:
            task_ids = self._extract_task_ids(target_count)
        else:
            task_ids = target_task_ids

        self.mc.log_debug('{} tasks will be inferred.'.format(len(task_ids)))

        in_queue = Queue()
        out_queue = Queue()

        procs = [Process(target=_dispatch_tasks, args=(task_ids, self.n_process, in_queue))]
        for i in range(self.n_process):
            procs.append(Process(target=_handle_tasks, args=(i, self.anno_path, self.bam_param,
                                                             self.param.mcmc_samples_count,
                                                             self.param.tune_samples_count,
                                                             self.theano_cache_dir,
                                                             in_queue, out_queue, self.mc)))

        for proc in procs:
            proc.start()

        self._store_tmp_result_records(out_queue)

        for proc in procs:
            proc.join()

        self._merge_results()
        self.mc.handle_progress('All tasks finished.')

        shutil.rmtree(self.tmp_dir)


def inference(args):
    """Inference."""
    anno_path = args['task']

    bam_paths = args['bam']
    read_len = args['read_len']

    paired_end = args['pe']
    insert_size_mean = args['insert_size_mean']
    insert_size_std = args['insert_size_std']

    out_dir = args['out_dir']

    n_process = args['process']

    mcmc_samples_count = 500
    tune_samples_count = 500
    if 'n_mcmc' in args:
        mcmc_samples_count = args['n_mcmc']
    if 'tune' in args:
        tune_samples_count = args['tune']

    target_count = args['count']

    mc = args['mc']

    bam_param = BAMParam(bam_paths=bam_paths, read_len=read_len, paired_end=paired_end,
                         insert_size_mean=insert_size_mean, insert_size_std=insert_size_std)
    param = InferenceParam(anno_path=anno_path, bam_param=bam_param,
                           mcmc_samples_count=mcmc_samples_count, tune_samples_count=tune_samples_count)

    inference_tool = InferenceTool(out_dir=out_dir, n_process=n_process, param=param, mc=mc)
    inference_tool.run(target_count)


def inference_resume(args):
    """Inference resume."""
    out_dir = args['out_dir']
    n_process = args['process']
    target_count = args['count']
    mc = args['mc']

    inference_tool = InferenceTool(out_dir=out_dir, n_process=n_process, param=None, mc=mc)
    inference_tool.run(target_count)
