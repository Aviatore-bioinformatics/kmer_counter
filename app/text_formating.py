COLOR_GREEN = '\033[92m'
COLOR_RED = '\033[31m'
COLOR_WHITE = '\033[0m'


def warning(msg):
    return f'{COLOR_RED}{msg}{COLOR_WHITE}'


def ok(msg):
    return f'{COLOR_GREEN}{msg}{COLOR_WHITE}'
