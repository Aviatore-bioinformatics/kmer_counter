COLOR_GREEN = '\033[92m'
COLOR_RED = '\033[91m'
COLOR_WHITE = '\033[0m'
COLOR_YELLOW = '\033[93m'


def red(msg):
    return f'{COLOR_RED}{msg}{COLOR_WHITE}'


def green(msg):
    return f'{COLOR_GREEN}{msg}{COLOR_WHITE}'


def print_warning(msg):
    print(f'[{COLOR_RED}{"WARNING"}{COLOR_WHITE}] - {msg}')


def print_info(msg, header=None):
    if header is None:
        print(f'[{COLOR_GREEN}{"INFO"}{COLOR_WHITE}] - {msg}', flush=True)
    else:
        print(f'[{COLOR_GREEN}{"INFO"}{COLOR_WHITE}][{COLOR_YELLOW}{header}{COLOR_WHITE}] - {msg}', flush=True)


def print_logo(msg):
    print("")
    print('#' * (len(msg) + 8))
    print(f"#{msg:^{len(msg) + 6}}#")
    print('#' * (len(msg) + 8))