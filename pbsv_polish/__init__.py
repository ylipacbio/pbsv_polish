VERSION = (0, 1, 0)


def get_version():
    return ".".join([str(i) for i in VERSION])


__version__ = get_version()


def f(a, b):
    return ' '.join([a, b])


PBSVP = 'pbsvp'
PBSVPUTIL = 'pbsvputil'

POLISH_ENTRY = f(PBSVP, 'polish')
COLLECT_ENTRY = f(PBSVP, 'collect')
TRIM_ENTRY = f(PBSVPUTIL, 'trim-lq')
SVDAGCON_ENTRY = f(PBSVPUTIL, 'svdagcon')
TRANSFORM_ENTRY = f(PBSVPUTIL, 'transform-coordinate')
