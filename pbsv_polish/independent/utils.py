#!/usr/bin/env python
from pbsv.independent.utils import _is_fmt, cmds_to_bash, execute, realpath, autofmt, is_fasta


def is_fastq(fn):
    """Return true if a file extension is fq or fastq"""
    return _is_fmt(fn, ["fq", "fastq"])


def get_movie_and_zmw_from_name(name):
    """Given a string of pacbio zmw name or read name, return movie and zmw
    ...doctest
        >>> get_movie_and_zmw_from_name('movie/1000/0_100')
        movie, 1000
    """
    try:
        fs = name.strip().split(' ')[0].split('/')
        movie, zmw = fs[0], fs[1]
        return movie, int(zmw)
    except ValueError:
        raise ValueError("Read %r is not a PacBio read." % name)


def get_zmws_from_iter(names_iter):
    """
    Return a list of zmw strings from an iterable of read names or zmws.
    ...doctest:
        >>> parse_zmws(['movie/0/0_1', 'movie2/2/10_11', 'movie/0/13_14'])
        ['movie2/2', 'movie/0']
    """
    zmws = set()
    for r in names_iter:
        movie, zmw = get_movie_and_zmw_from_name(r)
        zmws.add('{}/{}'.format(movie, zmw))
    return list(zmws)
