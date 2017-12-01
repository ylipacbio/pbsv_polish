import pytest

import pbsv_polish # test import
from pbsv_polish.cli import  pbsvp_main, pbsvputil_main


def g(f, args):
    with pytest.raises(SystemExit) as excinfo:
        f(args)
    assert 0 == excinfo.value.code

def test_pbsvp_main_align(capsys):
    g(pbsvp_main, ['main', 'polish', '--help'])
    g(pbsvp_main, ['main', 'collect', '-h'])

def test_pbsvputil_main_ngmlr(capsys):
    g(pbsvputil_main, ['main', 'trim-lq', '-h'])
    g(pbsvputil_main, ['main', 'svdagcon', '-h'])
    g(pbsvputil_main, ['main', 'transform-coordinate', '-h'])
