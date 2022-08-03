# MIT License
#
# Copyright (c) 2018-2021 Tskit Developers
# Copyright (c) 2015-2018 University of Oxford
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Python implementation of the simplify algorithm.
"""
import sys

import numpy as np
import portion

import tskit


def overlapping_segments(segments):
    """
    Returns an iterator over the (left, right, X) tuples describing the
    distinct overlapping segments in the specified set.
    """
    S = sorted(segments, key=lambda x: x.left)
    n = len(S)
    # Insert a sentinel at the end for convenience.
    S.append(Segment(sys.float_info.max, 0))
    right = S[0].left
    X = []
    j = 0
    while j < n:
        # Remove any elements of X with right <= left
        left = right
        X = [x for x in X if x.right > left]
        if len(X) == 0:
            left = S[j].left
        while j < n and S[j].left == left:
            X.append(S[j])
            j += 1
        j -= 1
        right = min(x.right for x in X)
        right = min(right, S[j + 1].left)
        yield left, right, X
        j += 1

    while len(X) > 0:
        left = right
        X = [x for x in X if x.right > left]
        if len(X) > 0:
            right = min(x.right for x in X)
            yield left, right, X


class Segment:
    """
    A class representing a single segment. Each segment has a left and right,
    denoting the loci over which it spans, a node and a next, giving the next
    in the chain.

    The node it records is the *output* node ID.
    """

    def __init__(self, left=None, right=None, node=None, next_segment=None):
        self.left = left
        self.right = right
        self.node = node
        self.next = next_segment

    def __str__(self):
        s = "({}-{}->{}:next={})".format(
            self.left, self.right, self.node, repr(self.next)
        )
        return s

    def __repr__(self):
        return repr((self.left, self.right, self.node))

    def __lt__(self, other):
        return (self.left, self.right, self.node) < (other.left, other.right, self.node)

class AncestorMap:
    """
    Simplifies a tree sequence to show relationships between
    samples and a designated set of ancestors.
    """

    def __init__(self, ts, sample, ancestors):
        self.ts = ts
        self.samples = set(sample)
        assert (self.samples).issubset(set(range(0, ts.num_nodes)))
        self.ancestors = set(ancestors)
        assert (self.ancestors).issubset(set(range(0, ts.num_nodes)))
        self.table = tskit.EdgeTable()
        self.sequence_length = ts.sequence_length
        self.A_head = [None for _ in range(ts.num_nodes)]
        self.A_tail = [None for _ in range(ts.num_nodes)]
        for sample_id in sample:
            self.add_ancestry(0, self.sequence_length, sample_id, sample_id)
        self.edge_buffer = {}

    def link_ancestors(self):
        if self.ts.num_edges > 0:
            all_edges = list(self.ts.edges())
            edges = all_edges[:1]
            for e in all_edges[1:]:
                if e.parent != edges[0].parent:
                    self.process_parent_edges(edges)
                    edges = []
                edges.append(e)
            self.process_parent_edges(edges)
        return self.table

    def process_parent_edges(self, edges):
        """
        Process all of the edges for a given parent.
        """
        assert len({e.parent for e in edges}) == 1
        parent = edges[0].parent
        S = []
        for edge in edges:
            x = self.A_head[edge.child]
            while x is not None:
                if x.right > edge.left and edge.right > x.left:
                    y = Segment(
                        max(x.left, edge.left), min(x.right, edge.right), x.node
                    )
                    S.append(y)
                x = x.next
        self.merge_labeled_ancestors(S, parent)
        self.check_state()

    def merge_labeled_ancestors(self, S, input_id):
        """
        All ancestry segments in S come together into a new parent.
        The new parent must be assigned and any overlapping segments coalesced.
        """
        is_sample = input_id in self.samples
        if is_sample:
            # Free up the existing ancestry mapping.
            x = self.A_tail[input_id]
            assert x.left == 0 and x.right == self.sequence_length
            self.A_tail[input_id] = None
            self.A_head[input_id] = None

        is_ancestor = input_id in self.ancestors
        prev_right = 0
        for left, right, X in overlapping_segments(S):
            if is_ancestor or is_sample:
                for x in X:
                    ancestry_node = x.node
                    self.record_edge(left, right, input_id, ancestry_node)
                self.add_ancestry(left, right, input_id, input_id)

                if is_sample and left != prev_right:
                    # Fill in any gaps in the ancestry for the sample.
                    self.add_ancestry(prev_right, left, input_id, input_id)

            else:
                for x in X:
                    ancestry_node = x.node
                    # Add sample ancestry for the currently-processed segment set.
                    self.add_ancestry(left, right, ancestry_node, input_id)
            prev_right = right

        if is_sample and prev_right != self.sequence_length:
            # If a trailing gap exists in the sample ancestry, fill it in.
            self.add_ancestry(prev_right, self.sequence_length, input_id, input_id)
        if input_id != -1:
            self.flush_edges()

    def record_edge(self, left, right, parent, child):
        """
        Adds an edge to the output list.
        """
        if child not in self.edge_buffer:
            self.edge_buffer[child] = [tskit.Edge(left, right, parent, child)]
        else:
            last = self.edge_buffer[child][-1]
            if last.right == left:
                last.right = right
            else:
                self.edge_buffer[child].append(tskit.Edge(left, right, parent, child))

    def add_ancestry(self, left, right, node, current_node):
        tail = self.A_tail[current_node]
        if tail is None:
            x = Segment(left, right, node)
            self.A_head[current_node] = x
            self.A_tail[current_node] = x
        else:
            if tail.right == left and tail.node == node:
                tail.right = right
            else:
                x = Segment(left, right, node)
                tail.next = x
                self.A_tail[current_node] = x

    def flush_edges(self):
        """
        Flush the edges to the output table after sorting and squashing
        any redundant records.
        """
        num_edges = 0
        for child in sorted(self.edge_buffer.keys()):
            for edge in self.edge_buffer[child]:
                self.table.append(edge)
                num_edges += 1
        self.edge_buffer.clear()
        return num_edges

    def check_state(self):

        num_nodes = len(self.A_head)
        for j in range(num_nodes):
            head = self.A_head[j]
            tail = self.A_tail[j]
            if head is None:
                assert tail is None
            else:
                x = head
                while x.next is not None:
                    x = x.next
                assert x == tail
                x = head.next
                while x is not None:
                    assert x.left < x.right
                    if x.next is not None:
                        if self.ancestors is None:
                            assert x.right <= x.next.left
                        # We should also not have any squashable segments.
                        if x.right == x.next.left:
                            assert x.node != x.next.node
                    x = x.next

    def print_state(self):
        print(".................")
        print("Ancestors: ")
        num_nodes = len(self.A_tail)
        for j in range(num_nodes):
            print("\t", j, "->", end="")
            x = self.A_head[j]
            while x is not None:
                print(f"({x.left}-{x.right}->{x.node})", end="")
                x = x.next
            print()
        print("Output:")
        print(self.table)
        self.check_state()


if __name__ == "__main__":
    # Simple CLI for running simplifier/ancestor mapping above.
    # class_to_implement = sys.argv[1]
    # assert class_to_implement == "Simplifier" or class_to_implement == "AncestorMap"
    ts = tskit.load(sys.argv[1])

    # elif class_to_implement == "AncestorMap":
    samples = ts.samples()

    census_time = sys.argv[2]
    # ancestors = ancestors.split(",")
    # ancestors = list(map(int, ancestors))
    ancestors = [u.id for u in ts.nodes() if u.time == census_time]

    s = AncestorMap(ts, samples, ancestors)
    tss = s.link_ancestors()
    # tables = tss.dump_tables()
    # print(tables.nodes)
    # print(tss)
