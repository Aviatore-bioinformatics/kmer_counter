COLOR_GREEN = '\033[92m'
COLOR_RED = '\033[31m'
COLOR_WHITE = '\033[0m'


def red(msg):
    return f'{COLOR_RED}{msg}{COLOR_WHITE}'


def green(msg):
    return f'{COLOR_GREEN}{msg}{COLOR_WHITE}'


def print_warning(msg):
    print(f'[{COLOR_RED}{"WARNING"}{COLOR_WHITE}] - {msg}')


def print_info(msg):
    print(f'[{COLOR_GREEN}{"INFO"}{COLOR_WHITE}] - {msg}')
