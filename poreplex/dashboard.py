#
# Copyright (c) 2018 Institute for Basic Science
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

from urwid.canvas import apply_text_layout
import urwid
import asyncio
import pandas as pd
from asyncio import CancelledError
from collections import Counter
from datetime import datetime
from functools import partial
import random
from random import choice
from weakref import proxy
from . import __version__


GROUP_NAME_LABELS = [
    (0, 'Barcode 1'),
    (1, 'Barcode 2'),
    (2, 'Barcode 3'),
    (3, 'Barcode 4'),
    (None, 'Undetermined'),
]
GROUP_NAME_UNBARCODED = 'All Reads'


class SelectGroupButton(urwid.Button):
    def __init__(self, caption):
        self.caption = caption
        cursor_position = len(caption) + 1
        super(SelectGroupButton, self).__init__("")
        self._w = urwid.AttrMap(urwid.SelectableIcon(
            ['', caption], cursor_position), None, 'selected')


class BidiJustifiedText(urwid.Text):

    def get_text(self, maxcol=None):
        text, attr = super(BidiJustifiedText, self).get_text()
        if maxcol is not None or self._cache_maxcol is not None:
            if maxcol is None:
                maxcol = self._cache_maxcol

            while maxcol < len(text):
                maxcol += maxcol
            fillersize = maxcol - len(text)
            if fillersize > 0:
                leftsize = sum(size for _, size in attr[:-1])
                text = text[:leftsize] + ' ' * fillersize + text[leftsize:]
                attr = attr[:]
                attr.insert(-1, (None, fillersize))
        return text, attr

    def _update_cache_translation(self,maxcol, ta):
        self._cache_maxcol = maxcol
        super(BidiJustifiedText, self)._update_cache_translation(maxcol, ta)


class ReadOnlyListBox(urwid.ListBox):
    _selectable = False


class ReadMappingStatistics:

    def __init__(self, groups):
        self.groups = groups
        self.total_counts = {}
        self.failed_counts = {}
        self.unmapped_counts = {}
        self.mapped_counts = {}
        self.update_queue = {}
        for group in list(groups) + ['_total']:
            self.total_counts[group] = \
            self.unmapped_counts[group] = \
            self.failed_counts[group] = 0
            self.mapped_counts[group] = Counter()
            self.update_queue[group] = []

    def update(self, group, mapped, unmapped, failed, update_total=True):
        self.total_counts[group] += len(mapped) + unmapped + failed
        if failed > 0:
            self.failed_counts[group] += failed
        if unmapped > 0:
            self.unmapped_counts[group] += unmapped
        if len(mapped) > 0:
            self.update_queue[group].extend(mapped)

        if update_total:
            self.update('_total', mapped, unmapped, failed, False)

    def most_common(self, group, num):
        counter = self.mapped_counts[group]
        queue = self.update_queue[group]
        if len(queue) > 0:
            counter.update(queue)
            del queue[:]

        return counter.most_common(num)

    def overall_stats(self, group):
        return (
            self.total_counts[group],
            self.failed_counts[group],
            self.unmapped_counts[group],
            self.total_counts[group] - self.failed_counts[group] - self.unmapped_counts[group],
        )

    def get_demux_rate(self):
        demuxed_count = sum(self.total_counts[grp] for grp in self.groups
                            if isinstance(grp, int))
        total_count = sum(self.total_counts[grp] for grp in self.groups)
        return demuxed_count / total_count if total_count > 0 else 0

    def get_mapped_rate(self):
        unmapped_count = sum(self.unmapped_counts[grp] for grp in self.groups)
        total_count = sum(self.total_counts[grp] for grp in self.groups)
        return (total_count - unmapped_count) / total_count if total_count > 0 else 0


class DashboardView:

    palette = [
        ('body',            'light gray',   'black'),
        ('header',          'dark magenta', 'white',        'bold'),
        ('footer',          'black',        'light gray'),
        ('overview label',  'dark cyan',    '',             'bold'),
        ('divider',         'dark gray',    ''),
        ('grpstatsname',    'brown',        ''),
        ('grpstatscount',   'light gray',   ''),
        ('ctgname',         'dark green',   ''),
        ('count',           'light gray',   ''),
        ('done1',           'white',        'dark gray'),
        ('complete1',       'black',        'light gray',   'bold'),
        ('done2',           'white',        'dark gray'),
        ('complete2',       'black',        'light gray',   'bold'),
        ('selected',        'black',        'dark cyan',    'standout'),
        ('quit button',     'dark red',     'light gray',   'bold'),
        ('arrow button',    'dark blue',    'light gray',   'bold'),
        ('modal frame',     'white',        'dark blue'),
        ('modal text',      'light gray',   'dark blue'),
    ]
    trstats_view_length = 50

    label_readstats_total = ('grpstatsname', 'Total')
    label_readstats_failed = ('grpstatsname', 'Failed')
    label_readstats_mapped = ('grpstatsname', 'Mapped')
    label_readstats_unmapped = ('grpstatsname', 'Non-mapped')

    label_overview_total = ('overview label', 'Total reads: ')
    label_overview_inqueue = ('overview label', 'In queue: ')
    label_overview_elapsed_time = ('overview label', 'Elapsed time: ')
    label_overview_active_tasks = ('overview label', 'Active tasks: ')

    overview_labels = {
        'demux_rate': 'Reads with identified barcodes',
        'mapped_rate': 'Mappable reads',
        'progress': 'Reads processed',
    }

    def __init__(self, session, readgroups, pbar_left_src, pbar_right_src,
                 analysis_delay=0, contigaliases={}):
        readgroups = self.adopt_readgroups(readgroups)
        self.session = proxy(session)
        self.urwid_loop = None
        self.readstats = ReadMappingStatistics([groupid for groupid, _ in readgroups])
        self.readgroups = readgroups
        self.contigaliases = contigaliases
        self.analysis_delay = analysis_delay

        self.readstats_left = self.readstats_right = None
        self.readstatsbox_left = None
        self.selected_group = readgroups[0][0]

        self.pbar_left = self.pbar_right = None
        self.pbar_left_value = self.pbar_right_value = 0
        self.pbar_left_source = pbar_left_src
        self.pbar_right_source = pbar_right_src
        self.pbar_left_getter = getattr(self, 'get_' + pbar_left_src)
        self.pbar_right_getter = getattr(self, 'get_' + pbar_right_src)

        self.overview_widget_total = \
        self.overview_widget_inqueue = \
        self.overview_widget_elapsed_time = \
        self.overview_widget_active_analysis = None

    def adopt_readgroups(self, readgroups):
        readgroups = [
            (groupid, name)
            for groupid, name in GROUP_NAME_LABELS
            if groupid in readgroups
        ]

        if not any(isinstance(groupid, int) for groupid, name in readgroups): # w/o barcoding
            readgroups = [(None, GROUP_NAME_UNBARCODED)]

        return readgroups

    def exit_on_q(self, key):
        if key in ('q', 'Q'):
            self.stop()

    def stop(self):
        self.urwid_loop.stop()
        self.session.dashboard = None
        self.session.stop('REQUEST')

    def build_widgets(self, will_align):
        mainbody = self.build_mainbody_widgets()
        header = self.build_header_widgets()
        footer = self.build_footer_widgets()

        frame = urwid.Frame(mainbody, header=header, footer=footer)
        frame = urwid.AttrWrap(frame, 'body')

        if self.analysis_delay > 0 and will_align:
            self.widget_under_modal = frame
            return self.add_wait_notice_widget(frame, self.analysis_delay)
        else:
            self.widget_under_modal = None
            return frame

    def add_wait_notice_widget(self, widget, delay_seconds):
        message = ("Waiting for the first update of results. In live mode, "
                   "an analysis starts on {} seconds later to compensate "
                   "many potential problems from the filesystem delays.".format(delay_seconds))
        lines = 5
        title = "Please wait for a while"

        waitbox = urwid.AttrWrap(urwid.LineBox(urwid.Filler(urwid.AttrWrap(
            urwid.Text(message), 'modal text')),
            title=title), 'modal frame')
        return urwid.Overlay(waitbox, widget, align='center', width=58, valign='middle',
                             height=lines)

    def build_read_stats_widget(self, title):
        # Make the upper part that shows total and unmapped reads
        genrow = lambda label: (
            BidiJustifiedText([('grpstatsname', label),
                               ('grpstatscount', '-')], wrap='any'))
        text_elements = {'total': genrow('Total'), 'unmapped': genrow('Non-mapped'),
                         'failed': genrow('Failed'), 'mapped': genrow('Mapped')}
        top_box = ReadOnlyListBox(
            urwid.SimpleListWalker([
                text_elements['total'],
                text_elements['mapped'],
                text_elements['unmapped'],
                text_elements['failed']]))

        text_elements['contigs'] = [BidiJustifiedText(['', ''], wrap='any')
                                    for i in range(self.trstats_view_length)]
        bottom_box = ReadOnlyListBox(urwid.SimpleListWalker(text_elements['contigs']))

        column_widget = (
            urwid.LineBox(
                urwid.Padding(
                    urwid.Pile([
                        ('fixed', 4, top_box),
                        ('flow', urwid.AttrWrap(urwid.Divider('-'), 'divider')),
                        bottom_box
                    ]),
                    left=1, right=1
                ),
                title=title
            )
        )
        return column_widget, text_elements

    def on_group_selection_changed(self, groupwalker):
        group_label = groupwalker[groupwalker.focus].caption
        self.readstatsbox_left.set_title(group_label)
        self.selected_group = self.readgroups[groupwalker.focus][0]
        self.update_group_read_stats(self.selected_group, self.readstats_left, clear=True)

    def build_mainbody_widgets(self):
        # Build the left-most group selector list box
        groupwalker = urwid.SimpleListWalker([
            SelectGroupButton(name) for _, name in self.readgroups])
        urwid.connect_signal(groupwalker, 'modified',
                             partial(self.on_group_selection_changed,
                                     groupwalker))

        group_selector = (
            urwid.Padding(
                urwid.Filler(
                    urwid.BoxAdapter(urwid.ListBox(groupwalker), 7),
                    'top', top=1
                ),
                left=1
            )
        )

        default_group_title = self.readgroups[0][1]
        left_widget, left_elements = self.build_read_stats_widget(default_group_title)
        right_widget, right_elements = self.build_read_stats_widget('All Reads')

        stats_columns = urwid.Columns([
            ('fixed', 16, group_selector),
            left_widget, right_widget,
        ])

        self.readstats_left = left_elements
        self.readstats_right = right_elements

        self.readstatsbox_left = left_widget

        return urwid.Filler(stats_columns, 'top', height=('relative', 100))

    def build_header_widgets(self):
        # Top title bar
        titlerow = urwid.Text(' Poreplex {}'.format(__version__), align='left')
        titlerow = urwid.AttrWrap(titlerow, 'header')

        # Left status column
        self.pbar_left = urwid.ProgressBar('done1', 'complete1', 0, 100)

        self.overview_widget_total = urwid.Text([self.label_overview_total, '0'])
        self.overview_widget_inqueue = urwid.Text([self.label_overview_inqueue, '0'])

        overview_left = urwid.Pile([
            self.overview_widget_total,
            self.overview_widget_inqueue,
            urwid.Text([('overview label',
                        self.overview_labels[self.pbar_left_source] + ':')]),
            self.pbar_left,
        ])

        # Right status column
        self.pbar_right = urwid.ProgressBar('done2', 'complete2', 0, 100)

        self.overview_widget_elapsed_time = urwid.Text(
            [self.label_overview_elapsed_time, '0:00:00'])
        self.overview_widget_active_analysis = urwid.Text(
            [self.label_overview_active_tasks, '0'])

        overview_right = urwid.Pile([
            self.overview_widget_elapsed_time,
            self.overview_widget_active_analysis,
            urwid.Text([('overview label',
                        self.overview_labels[self.pbar_right_source] + ':')]),
            self.pbar_right,
        ])

        # Container for the status columns
        overview = urwid.Columns([
            urwid.Padding(overview_left, 'center', left=1, right=1),
            urwid.Padding(overview_right, 'center', left=1, right=1),
        ])
        overview = urwid.Filler(overview, 'middle', height='pack', top=1, bottom=1)
        overview = urwid.BoxAdapter(overview, 6)

        return urwid.Pile([titlerow, overview])

    def build_footer_widgets(self):
        instruction_line = urwid.Text([
            'Press (', ('arrow button', 'Up'), ' or ',
            ('arrow button', 'Down'), ') to select a group to view. ',
            'Press (', ('quit button', 'Q'), ') to stop and quit.'
        ])
        return urwid.AttrWrap(instruction_line, 'footer')

    def start(self, loop, will_align):
        widgets = self.build_widgets(will_align)

        uaioloop = urwid.AsyncioEventLoop(loop=loop)
        self.urwid_loop = urwid.MainLoop(
            widgets, event_loop=uaioloop,
            palette=self.palette,
            unhandled_input=self.exit_on_q
        )

        self.start_update_hooks(loop)

        self.urwid_loop.start()
        #uaioloop.run()

    def start_update_hooks(self, loop):
        loop.create_task(self.update_elapsed_time())
        loop.create_task(self.update_read_count_stats())
        loop.create_task(self.update_overview_status())

    async def update_elapsed_time(self):
        label = self.label_overview_elapsed_time
        widget = self.overview_widget_elapsed_time
        start_time = datetime.now()

        try:
            while True:
                timediff = str(datetime.now() - start_time).split('.')[0]
                widget.set_text([label, timediff])
                await asyncio.sleep(1)
        except CancelledError:
            pass

    async def update_read_count_stats(self):
        try:
            while True:
                self.update_group_read_stats('_total', self.readstats_right)
                self.update_group_read_stats(self.selected_group, self.readstats_left)

                await asyncio.sleep(2)
        except CancelledError:
            pass

    def get_demux_rate(self):
        return int(self.readstats.get_demux_rate() * 100)

    def get_mapped_rate(self):
        return int(self.readstats.get_mapped_rate() * 100)

    def get_progress(self):
        if self.session.reads_found > 0:
            return int(self.session.reads_processed / self.session.reads_found * 100)
        else:
            return 0

    async def update_overview_status(self):
        try:
            while True:
                # Left column stats
                self.overview_widget_total.set_text([
                    self.label_overview_total,
                    format(self.session.reads_found, ',d')])
                self.overview_widget_inqueue.set_text([
                    self.label_overview_inqueue,
                    format(self.session.reads_queued, ',d')])

                # Left progress bar
                leftnew = self.pbar_left_getter()
                if leftnew != self.pbar_left_value:
                    self.pbar_left.set_completion(leftnew)
                    self.pbar_left_value = leftnew

                # Right column stats
                self.overview_widget_active_analysis.set_text([
                    self.label_overview_active_tasks,
                    str(self.session.active_batches)])

                # Right progress bar
                rightnew = self.pbar_right_getter()
                if rightnew != self.pbar_right_value:
                    self.pbar_right.set_completion(rightnew)
                    self.pbar_right_value = rightnew

                await asyncio.sleep(1)
        except CancelledError:
            pass

    def update_group_read_stats(self, groupid, widgets, clear=False):
        toplist = self.readstats.most_common(groupid, self.trstats_view_length)
        total, failed, unmapped, mapped = self.readstats.overall_stats(groupid)

        widgets['total'].set_text([self.label_readstats_total,
                                   ('grpstatscount', format(total, ',d'))])
        widgets['mapped'].set_text([self.label_readstats_mapped,
                                   ('grpstatscount', format(mapped, ',d'))])
        widgets['failed'].set_text([self.label_readstats_failed,
                                    ('grpstatscount', format(failed, ',d'))])
        widgets['unmapped'].set_text([self.label_readstats_unmapped,
                                      ('grpstatscount', format(unmapped, ',d'))])

        aliases = self.contigaliases
        for (trname, count), widget in zip(toplist, widgets['contigs']):
            widget.set_text([('ctgname', aliases.get(trname, trname)),
                             ('count', format(count, ',d'))])

        if clear and len(toplist) < len(widgets['contigs']):
            for i in range(len(toplist), len(widgets['contigs'])):
                widgets['contigs'][i].set_text(['', ''])

    def hide_modal(self):
        self.urwid_loop.widget = self.widget_under_modal
        self.widget_under_modal = None

    def feed_mapped(self, mapresult):
        for rgid, _ in self.readgroups:
            self.readstats.update(rgid, mapresult['mapped'][rgid],
                                  mapresult['unmapped'][rgid],
                                  mapresult['failed'][rgid])

        if self.widget_under_modal is not None:
            self.hide_modal()

def load_aliases(filename):
    return (
        pd.read_table(filename, names=['contig', 'alias'], index_col=0)
            ['alias'].to_dict()
    )

